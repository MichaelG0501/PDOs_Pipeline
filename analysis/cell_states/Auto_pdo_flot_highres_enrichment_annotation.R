####################
# Auto_pdo_flot_highres_enrichment_annotation.R
#
# Enrichment annotation for matched-FLOT high-resolution PDO metaprograms.
# Consumes retained MP genes from Auto_pdo_flot_matched_highres_mp_trend_filter.R
# and plots all reference enrichment families as a multi-page PDF plus PNG pages.
#
# Env: dmtcp
####################

args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(grid)
})

####################
# setup
####################
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_pdo_flot_highres_metaprogram_trends"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

type_levels <- c("increase", "decrease")
patient_order <- c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
sample_order <- as.vector(rbind(
  paste0(patient_order, "_Untreated_PDO"),
  paste0(patient_order, "_Treated_PDO")
))
max_mps_per_enrichment_page <- 18
enrichment_cap <- 7

####################
# helpers
####################
safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

chunk_vector <- function(x, n) {
  if (length(x) == 0) return(list())
  split(x, ceiling(seq_along(x) / n))
}

write_mp_gene_table <- function(mp_genes, path) {
  if (length(mp_genes) == 0) {
    gene_table <- data.frame(MP = character(), rank = integer(), gene = character(), stringsAsFactors = FALSE)
  } else {
    gene_table <- do.call(rbind, lapply(names(mp_genes), function(mp) {
      data.frame(MP = mp, rank = seq_along(mp_genes[[mp]]), gene = mp_genes[[mp]], stringsAsFactors = FALSE)
    }))
  }
  write.csv(gene_table, path, row.names = FALSE)
}

order_mps_by_metrics <- function(mp_vec, trend_summary) {
  if (length(mp_vec) == 0) return(character())
  trend_summary |>
    dplyr::filter(MP %in% mp_vec) |>
    dplyr::mutate(
      numberPrograms = suppressWarnings(as.integer(numberPrograms)),
      numberPrograms = tidyr::replace_na(numberPrograms, -Inf),
      silhouette = tidyr::replace_na(silhouette, -Inf),
      pair_support_n = tidyr::replace_na(pair_support_n, -Inf),
      trend_p_value = tidyr::replace_na(trend_p_value, Inf),
      mean_delta_min_abs = tidyr::replace_na(mean_delta_min_abs, -Inf)
    ) |>
    dplyr::arrange(trend_p_value, dplyr::desc(pair_support_n), dplyr::desc(mean_delta_min_abs), dplyr::desc(numberPrograms), dplyr::desc(silhouette), MP) |>
    dplyr::pull(MP)
}

order_retained_by_type <- function(trend_summary) {
  unlist(lapply(type_levels, function(direction_name) {
    order_mps_by_metrics(trend_summary$MP[trend_summary$retained & trend_summary$treatment_direction == direction_name], trend_summary)
  }), use.names = FALSE)
}

make_display_labels <- function(trend_summary, nMP) {
  label_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_top_3CA_noncellcycle_nMP", nMP, ".csv"))
  label_table <- if (file.exists(label_path)) {
    read.csv(label_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    data.frame(MP = trend_summary$MP, top_3ca_noncc = NA_character_, stringsAsFactors = FALSE)
  }
  label_df <- trend_summary |>
    dplyr::left_join(label_table, by = "MP") |>
    dplyr::mutate(
      top_3ca_noncc = ifelse(is.na(top_3ca_noncc) | top_3ca_noncc == "", "3CA:no_nonCC_hit", top_3ca_noncc),
      single_programme_label = ifelse(numberPrograms == 1, "\n[1 programme]", ""),
      display_label = paste0(MP, "\n", top_3ca_noncc, single_programme_label)
    )
  setNames(label_df$display_label, label_df$MP)
}

load_enrichment_references <- function() {
  hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
  hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

  mp_csv <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
  if (file.exists(mp_csv)) {
    mp_list <- read.csv(mp_csv, check.names = FALSE, stringsAsFactors = FALSE)
    mp_list <- as.list(mp_list)
    mp_list <- lapply(mp_list, function(x) unique(x[x != "" & !is.na(x)]))
    mp_term2gene <- data.frame(
      term = rep(names(mp_list), lengths(mp_list)),
      gene = unlist(mp_list, use.names = FALSE),
      row.names = NULL,
      stringsAsFactors = FALSE
    )
    mp_term2gene$term <- sub("^MP", "3CA_mp", mp_term2gene$term)
    mp_term2name <- data.frame(term = unique(mp_term2gene$term), name = unique(mp_term2gene$term), stringsAsFactors = FALSE)
  } else {
    warning("3CA MP reference not found: ", mp_csv)
    mp_term2gene <- data.frame(term = character(), gene = character(), stringsAsFactors = FALSE)
    mp_term2name <- data.frame(term = character(), name = character(), stringsAsFactors = FALSE)
  }

  individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
  if (dir.exists(individual_dir)) {
    custom_files <- list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE)
    custom_refs <- lapply(custom_files, readRDS)
    names(custom_refs) <- sub("\\.rds$", "", sub(".*enrich_dev_", "", basename(custom_files)))
  } else {
    warning("Developmental reference directory not found: ", individual_dir)
    custom_refs <- list()
  }

  list(
    hallmark_term2gene = hallmark_term2gene,
    hallmark_term2name = hallmark_term2name,
    mp_term2gene = mp_term2gene,
    mp_term2name = mp_term2name,
    custom_refs = custom_refs
  )
}

extract_enrichment_rows <- function(cluster_enrich, element, include_non_significant = FALSE) {
  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)
    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    r <- r[is.finite(r$p.adjust) & r$p.adjust > 0, , drop = FALSE]
    if (!include_non_significant) {
      r <- r[r$p.adjust < 0.05, , drop = FALSE]
    }
    if (nrow(r) == 0) return(NULL)
    term <- if ("Description" %in% colnames(r)) r$Description else r$ID
    data.frame(
      Program = prog,
      Term = term,
      padj = r$p.adjust,
      Overlap = r$GeneRatio,
      Count = if ("Count" %in% colnames(r)) r$Count else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(df_list)
}

choose_terms_for_chunk <- function(df, mp_chunk, top_per_program = 8, top_n = 80) {
  terms <- df |>
    dplyr::filter(Program %in% mp_chunk, padj < 0.05) |>
    dplyr::arrange(match(Program, mp_chunk), padj, dplyr::desc(Count)) |>
    dplyr::group_by(Program) |>
    dplyr::slice_head(n = top_per_program) |>
    dplyr::ungroup() |>
    dplyr::distinct(Term) |>
    dplyr::pull(Term)

  if (length(terms) > top_n) {
    terms <- df |>
      dplyr::filter(Program %in% mp_chunk, Term %in% terms) |>
      dplyr::group_by(Term) |>
      dplyr::summarise(min_p = min(padj, na.rm = TRUE), max_count = max(Count, na.rm = TRUE), .groups = "drop") |>
      dplyr::arrange(min_p, dplyr::desc(max_count)) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::pull(Term)
  }
  terms
}

order_rows_by_significance <- function(mat, term_padj, mp_chunk) {
  if (nrow(mat) <= 1) return(rownames(mat))
  best_idx <- max.col(mat, ties.method = "first")
  best_program <- colnames(mat)[best_idx]
  row_order <- unlist(lapply(mp_chunk[mp_chunk %in% colnames(mat)], function(mp) {
    rows <- rownames(mat)[best_program == mp]
    if (length(rows) == 0) return(character())
    if (length(rows) > 2) {
      rows[hclust(dist(mat[rows, , drop = FALSE]), method = "ward.D2")$order]
    } else {
      rows[order(term_padj[rows], decreasing = FALSE, na.last = TRUE)]
    }
  }), use.names = FALSE)
  c(row_order, setdiff(rownames(mat), row_order))
}

draw_enrichment_heatmap <- function(mat, text_mat, trend_summary, label_map, element, page_label,
                                    png_path = NULL, draw_to_current_device = TRUE) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(FALSE)
  cols <- grDevices::colorRampPalette(c("white", "#fee0d2", "#fc9272", "#de2d26", "#67000d"))(100)
  breaks <- seq(0, enrichment_cap, length.out = length(cols) + 1)

  col_anno <- trend_summary |>
    dplyr::filter(MP %in% colnames(mat)) |>
    dplyr::mutate(Pair_support = paste0(pair_support_n, "/4")) |>
    dplyr::select(MP, Direction = treatment_direction, Pair_support) |>
    as.data.frame()
  row.names(col_anno) <- col_anno$MP
  col_anno <- col_anno[colnames(mat), c("Direction", "Pair_support"), drop = FALSE]
  col_gaps <- which(col_anno$Direction[-nrow(col_anno)] != col_anno$Direction[-1])
  anno_colors <- list(
    Direction = c(increase = "#4a4a4a", decrease = "#a0a0a0"),
    Pair_support = c("3/4" = "#c2855a", "4/4" = "#2b6a8e")
  )
  best_program <- colnames(mat)[max.col(mat, ties.method = "first")]
  row_gaps <- which(best_program[-length(best_program)] != best_program[-1])
  col_labels <- label_map[colnames(mat)]

  make_plot <- function(silent) {
    pheatmap::pheatmap(
      mat,
      display_numbers = text_mat,
      number_color = "black",
      fontsize_number = 7,
      labels_col = col_labels,
      color = cols,
      breaks = breaks,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      gaps_row = row_gaps,
      gaps_col = col_gaps,
      annotation_col = col_anno,
      annotation_colors = anno_colors,
      border_color = NA,
      show_colnames = TRUE,
      angle_col = 45,
      fontsize = 18,
      fontsize_row = 8,
      fontsize_col = 8,
      main = paste0(element, " Enrichment (-log10 adjusted p), ", page_label),
      silent = silent
    )
  }

  if (draw_to_current_device) {
    grid::grid.newpage()
    p_obj <- make_plot(silent = TRUE)
    grid::grid.draw(p_obj$gtable)
  }

  if (!is.null(png_path)) {
    png(png_path, width = 4200, height = max(2400, 85 * nrow(mat) + 900), res = 300)
    make_plot(silent = FALSE)
    dev.off()
  }

  TRUE
}

run_enrichment_heatmaps <- function(cluster_enrich, selected_mp_genes, trend_summary, label_map, element,
                                    output_prefix, top_per_program = 8, top_n = 80,
                                    write_png = TRUE, draw_pdf = TRUE) {
  df <- extract_enrichment_rows(cluster_enrich, element, include_non_significant = FALSE)
  if (is.null(df) || nrow(df) == 0) {
    message("No significant enrichment results to plot for: ", element)
    return(FALSE)
  }

  ordered_mps <- names(selected_mp_genes)
  chunks <- chunk_vector(ordered_mps, max_mps_per_enrichment_page)
  plotted <- logical(length(chunks))

  for (i in seq_along(chunks)) {
    mp_chunk <- chunks[[i]]
    terms_use <- choose_terms_for_chunk(df, mp_chunk, top_per_program = top_per_program, top_n = top_n)
    if (length(terms_use) == 0) next

    full_grid <- expand.grid(Term = terms_use, Program = mp_chunk, stringsAsFactors = FALSE)
    final_df <- full_grid |>
      dplyr::left_join(df, by = c("Term", "Program")) |>
      dplyr::mutate(
        score = tidyr::replace_na(pmin(-log10(padj), enrichment_cap), 0),
        display_text = tidyr::replace_na(Overlap, "")
      )

    mat <- final_df |>
      dplyr::select(Term, Program, score) |>
      tidyr::pivot_wider(names_from = Program, values_from = score) |>
      as.data.frame()
    row.names(mat) <- mat$Term
    mat <- as.matrix(dplyr::select(mat, -Term))
    mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

    text_mat <- final_df |>
      dplyr::select(Term, Program, display_text) |>
      tidyr::pivot_wider(names_from = Program, values_from = display_text) |>
      as.data.frame()
    row.names(text_mat) <- text_mat$Term
    text_mat <- as.matrix(dplyr::select(text_mat, -Term))

    term_padj <- df |>
      dplyr::filter(Term %in% rownames(mat), Program %in% colnames(mat)) |>
      dplyr::group_by(Term) |>
      dplyr::summarise(min_p = min(padj, na.rm = TRUE), .groups = "drop")
    term_padj <- setNames(term_padj$min_p, term_padj$Term)
    row_order <- order_rows_by_significance(mat, term_padj, mp_chunk)
    mat <- mat[row_order, , drop = FALSE]
    text_mat <- text_mat[row_order, , drop = FALSE]

    png_path <- if (write_png) {
      file.path(out_dir, paste0(output_prefix, "_", safe_name(element), "_group", i, ".png"))
    } else {
      NULL
    }
    plotted[i] <- draw_enrichment_heatmap(
      mat,
      text_mat,
      trend_summary,
      label_map,
      element,
      paste0("MP group ", i, "/", length(chunks)),
      png_path = png_path,
      draw_to_current_device = draw_pdf
    )
  }

  any(plotted)
}

safe_enrich <- function(expr) {
  tryCatch(expr, error = function(e) {
    warning(conditionMessage(e))
    NULL
  })
}

####################
# load retained MPs
####################
config_files <- list.files(out_dir, pattern = "^Auto_pdo_flot_highres_nMP[0-9]+_config\\.csv$", full.names = TRUE)
if (length(config_files) == 0) {
  stop("No high-resolution config found in ", out_dir, ". Run Auto_pdo_flot_matched_highres_mp_trend_filter.R first.")
}
config <- read.csv(config_files[order(file.info(config_files)$mtime, decreasing = TRUE)[1]], check.names = FALSE)
nMP <- as.integer(config$nMP[1])

selected_genes_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_selected_mp_genes_nMP", nMP, ".rds"))
trend_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_trend_summary_nMP", nMP, ".csv"))
if (!file.exists(selected_genes_path) || !file.exists(trend_path)) {
  stop("Missing selected MP outputs. Run Auto_pdo_flot_matched_highres_mp_trend_filter.R first.")
}

selected_mp_genes <- readRDS(selected_genes_path)
trend_summary <- read.csv(trend_path, check.names = FALSE, stringsAsFactors = FALSE)
trend_summary$retained <- trend_summary$retained == TRUE | trend_summary$retained == "TRUE"
retained_order <- order_retained_by_type(trend_summary)
selected_mp_genes <- selected_mp_genes[retained_order]
selected_mp_genes <- selected_mp_genes[lengths(selected_mp_genes) > 0]

if (length(selected_mp_genes) == 0) {
  writeLines("No retained MPs; enrichment annotation was not run.", file.path(out_dir, paste0("Auto_pdo_flot_highres_enrichment_skipped_nMP", nMP, ".txt")))
  message("No retained MPs; skipping enrichment.")
  quit(save = "no")
}

write_mp_gene_table(selected_mp_genes, file.path(out_dir, paste0("Auto_pdo_flot_highres_enrichment_input_mp_genes_nMP", nMP, ".csv")))
label_map <- make_display_labels(trend_summary, nMP)
refs <- load_enrichment_references()

####################
# run enrichment
####################
cluster_enrich_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_cluster_enrich_nMP", nMP, ".rds"))
if (file.exists(cluster_enrich_path) && !force) {
  message("Loading existing enrichment object: ", cluster_enrich_path)
  cluster_enrich <- readRDS(cluster_enrich_path)
} else {
  message("Running enrichment for ", length(selected_mp_genes), " retained MPs.")
  cluster_enrich <- lapply(names(selected_mp_genes), function(mp_name) {
    genes <- unique(selected_mp_genes[[mp_name]])
    message("Processing enrichment for ", mp_name)
    res_GO <- safe_enrich(clusterProfiler::enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      qvalueCutoff = 0.05,
      readable = TRUE
    ))
    res_H <- safe_enrich(clusterProfiler::enricher(
      gene = genes,
      TERM2GENE = refs$hallmark_term2gene,
      TERM2NAME = refs$hallmark_term2name,
      qvalueCutoff = 0.05
    ))
    res_M <- if (nrow(refs$mp_term2gene) > 0) {
      safe_enrich(clusterProfiler::enricher(
        gene = genes,
        TERM2GENE = refs$mp_term2gene,
        TERM2NAME = refs$mp_term2name,
        qvalueCutoff = 0.05
      ))
    } else {
      NULL
    }
    res_custom_list <- lapply(names(refs$custom_refs), function(ref_name) {
      message("  -> Running custom enrichment: ", ref_name)
      safe_enrich(clusterProfiler::enricher(
        gene = genes,
        TERM2GENE = refs$custom_refs[[ref_name]]$TERM2GENE,
        TERM2NAME = refs$custom_refs[[ref_name]]$TERM2NAME,
        pAdjustMethod = "BH",
        qvalueCutoff = 0.05
      ))
    })
    names(res_custom_list) <- names(refs$custom_refs)
    c(list(rep_prog = mp_name, genes = genes, GO = res_GO, Hallmark = res_H, MPs_3CA = res_M), res_custom_list)
  })
  names(cluster_enrich) <- names(selected_mp_genes)
  saveRDS(cluster_enrich, cluster_enrich_path, compress = FALSE)
}

elements <- c(
  "MPs_3CA",
  "Hallmark",
  "GO",
  intersect(
    c(
      "Early_Embryogenesis",
      "Organogenesis_major",
      "Organogenesis_sub",
      "Adult_Epithelium",
      "Barretts_Oesophagus",
      "Normal_Development_long",
      "Normal_Development_short"
    ),
    names(refs$custom_refs)
  )
)

flat_rows <- lapply(elements, function(element) {
  df <- extract_enrichment_rows(cluster_enrich, element, include_non_significant = TRUE)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df$reference <- element
  df
})
flat_enrich <- dplyr::bind_rows(flat_rows)
if (nrow(flat_enrich) > 0) {
  flat_enrich <- flat_enrich |>
    dplyr::left_join(trend_summary[, c("MP", "treatment_direction", "pair_support_n", "trend_p_value", "trend_p_adj")], by = c("Program" = "MP")) |>
    dplyr::arrange(reference, match(Program, names(selected_mp_genes)), padj)
} else {
  flat_enrich <- data.frame(
    Program = character(),
    Term = character(),
    padj = numeric(),
    Overlap = character(),
    Count = integer(),
    reference = character(),
    treatment_direction = character(),
    pair_support_n = integer(),
    trend_p_value = numeric(),
    trend_p_adj = numeric(),
    stringsAsFactors = FALSE
  )
}
write.csv(flat_enrich, file.path(out_dir, paste0("Auto_pdo_flot_highres_enrichment_results_nMP", nMP, ".csv")), row.names = FALSE)

####################
# plot enrichment
####################
output_prefix <- paste0("Auto_pdo_flot_highres_enrich_nMP", nMP)
pdf_path <- file.path(out_dir, paste0("Auto_pdo_flot_highres_enrichment_annotation_nMP", nMP, ".pdf"))
pdf(pdf_path, width = 20, height = 14, useDingbats = FALSE)
plotted <- logical(0)
plotted <- c(plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, "MPs_3CA", output_prefix, top_per_program = 8, top_n = 80, write_png = FALSE, draw_pdf = TRUE))
plotted <- c(plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, "Hallmark", output_prefix, top_per_program = 8, top_n = 80, write_png = FALSE, draw_pdf = TRUE))
plotted <- c(plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, "GO", output_prefix, top_per_program = 6, top_n = 60, write_png = FALSE, draw_pdf = TRUE))

for (element in setdiff(elements, c("MPs_3CA", "Hallmark", "GO"))) {
  plotted <- c(plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, element, output_prefix, top_per_program = 8, top_n = 80, write_png = FALSE, draw_pdf = TRUE))
}
dev.off()

png_plotted <- logical(0)
png_plotted <- c(png_plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, "MPs_3CA", output_prefix, top_per_program = 8, top_n = 80, write_png = TRUE, draw_pdf = FALSE))
png_plotted <- c(png_plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, "Hallmark", output_prefix, top_per_program = 8, top_n = 80, write_png = TRUE, draw_pdf = FALSE))
png_plotted <- c(png_plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, "GO", output_prefix, top_per_program = 6, top_n = 60, write_png = TRUE, draw_pdf = FALSE))
for (element in setdiff(elements, c("MPs_3CA", "Hallmark", "GO"))) {
  png_plotted <- c(png_plotted, run_enrichment_heatmaps(cluster_enrich, selected_mp_genes, trend_summary, label_map, element, output_prefix, top_per_program = 8, top_n = 80, write_png = TRUE, draw_pdf = FALSE))
}

if (!any(plotted)) {
  writeLines(
    "Enrichment ran, but no reference heatmaps had plottable significant content.",
    file.path(out_dir, paste0("Auto_pdo_flot_highres_enrichment_no_plottable_heatmaps_nMP", nMP, ".txt"))
  )
}

message("Enrichment outputs written to: ", out_dir)
message("Auto_pdo_flot_highres_enrichment_annotation.R completed successfully.")
