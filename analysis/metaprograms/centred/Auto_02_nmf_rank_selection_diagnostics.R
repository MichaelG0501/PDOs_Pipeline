####################
# Analysis registry:
#   Status: active
#   Script: analysis/metaprograms/centred/Auto_02_nmf_rank_selection_diagnostics.R
#   Methodology: analysis/methodology/metaprograms/centred/Auto_centred_metaprogram_refinement_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Computes silhouette + WSS diagnostics for centred NMF metaprograms across
#     a range of nMP values. Identifies optimal nMP via kneedle inflection-point
#     algorithm. Runs initial enrichment annotation on the optimal nMP result.
#   Inputs:
#     - PDOs_outs/centred_mp_refinement/geneNMF_metaprograms_nMP_{k}.rds
#   Outputs:
#     - live: PDOs_outs/centred_mp_refinement/optimal_nMP.rds
#     - live: PDOs_outs/centred_mp_refinement/figures/rank_selection_diagnostics_centred.pdf
#     - live: PDOs_outs/centred_mp_refinement/figures/centred_initial_enrichment_anno.pdf
#     - live: PDOs_outs/centred_mp_refinement/figures/Auto_centred_nMP_{optimal}_custom_heatmap_by_batch.pdf
#     - live: PDOs_outs/centred_mp_refinement/tables/Auto_centred_nMP_{optimal}_program_batch_annotations.csv
#   Downstream use:
#     - optimal_nMP.rds is an input to centred refinement steps 03 and 04.
#     - diagnostics, enrichment, and customized heatmap are terminal QC figures.
#   Conda env: dmtcp
####################

library(cluster)
library(ggplot2)
library(patchwork)

# === Paths ===
live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
outdir <- file.path(live_base, "centred_mp_refinement")
fig_dir <- file.path(outdir, "figures")
table_dir <- file.path(outdir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# === Compute metrics across nMP range ===
k_vals <- 4:25
avg_sil_widths <- numeric(length(k_vals))
wss_vals <- numeric(length(k_vals))

for (i in seq_along(k_vals)) {
  k <- k_vals[i]
  rds_path <- file.path(outdir, paste0("geneNMF_metaprograms_nMP_", k, ".rds"))

  if (file.exists(rds_path)) {
    mp_res <- readRDS(rds_path)

    # Pull distance matrix directly used during clustering
    dist_mat <- as.dist(1 - mp_res$programs.similarity)
    cluster_assignments <- cutree(mp_res$programs.tree, k = k)

    # 1. Silhouette score
    sil <- silhouette(cluster_assignments, dist = dist_mat)
    avg_sil_widths[i] <- summary(sil)$avg.width

    # 2. Within-Cluster Sum of Squares (WSS) based on Cosine Distance
    wss_k <- 0
    dist_m <- as.matrix(dist_mat)
    for (clust_id in unique(cluster_assignments)) {
      idx <- which(cluster_assignments == clust_id)
      if (length(idx) > 1) {
        cluster_dist <- dist_m[idx, idx]
        wss_k <- wss_k + sum(cluster_dist^2) / (2 * length(idx))
      }
    }
    wss_vals[i] <- wss_k
  } else {
    avg_sil_widths[i] <- NA
    wss_vals[i] <- NA
  }
}

df_metrics <- data.frame(nMP = k_vals, Silhouette = avg_sil_widths, WSS = wss_vals)
print(df_metrics)

# === Kneedle algorithm ===
find_knee <- function(x, y) {
  valid <- !is.na(y)
  x <- x[valid]; y <- y[valid]
  if (length(x) < 3) return(x[which.max(y)])
  x_norm <- (x - min(x)) / (max(x) - min(x))
  y_norm <- (y - min(y)) / (max(y) - min(y))
  x1 <- x_norm[1]; y1 <- y_norm[1]
  x2 <- x_norm[length(x_norm)]; y2 <- y_norm[length(y_norm)]
  dists <- abs((y2 - y1) * x_norm - (x2 - x1) * y_norm + x2 * y1 - y2 * x1) /
           sqrt((y2 - y1)^2 + (x2 - x1)^2)
  return(x[which.max(dists)])
}

sil_knee <- find_knee(df_metrics$nMP, df_metrics$Silhouette)
wss_knee <- find_knee(df_metrics$nMP, df_metrics$WSS)
optimal_nMP <- sil_knee

message(paste0("Silhouette inflection point: nMP = ", sil_knee))
message(paste0("WSS elbow point: nMP = ", wss_knee))
message(paste0("Selected optimal nMP: ", optimal_nMP))

saveRDS(optimal_nMP, file.path(outdir, "optimal_nMP.rds"))

# === Diagnostic plots ===
p1 <- ggplot(df_metrics, aes(x = nMP, y = Silhouette)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 3) +
  geom_vline(xintercept = sil_knee, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = sil_knee + 0.5, y = max(df_metrics$Silhouette, na.rm = TRUE),
           label = paste0("Inflection: ", sil_knee), hjust = 0, color = "red", size = 3.5) +
  theme_minimal() +
  scale_x_continuous(breaks = k_vals) +
  labs(title = "Silhouette Analysis (Centred PDOs)",
       x = "Number of MetaPrograms (nMP)",
       y = "Average Silhouette Width") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(df_metrics, aes(x = nMP, y = WSS)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "darkred", size = 3) +
  geom_vline(xintercept = wss_knee, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = wss_knee + 0.5, y = max(df_metrics$WSS, na.rm = TRUE) * 0.95,
           label = paste0("Elbow: ", wss_knee), hjust = 0, color = "red", size = 3.5) +
  theme_minimal() +
  scale_x_continuous(breaks = k_vals) +
  labs(title = "Elbow Method (WSS) (Centred PDOs)",
       x = "Number of MetaPrograms (nMP)",
       y = "Total Within Sum of Squares") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf(file.path(fig_dir, "rank_selection_diagnostics_centred.pdf"), width = 12, height = 6)
print(p1 + p2)
dev.off()

# ============================================================================
# Initial enrichment annotation on optimal nMP
# ============================================================================
cat("Running initial enrichment annotation...\n")
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(dplyr)
library(tidyr)
library(pheatmap)

geneNMF.metaprograms <- readRDS(
  file.path(outdir, paste0("geneNMF_metaprograms_nMP_", optimal_nMP, ".rds"))
)

run_enrichment_and_plot <- function(mp_list, valid_cluster_ids, mp_tree_order, out_pdf, cols_palette) {
  hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
  hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
  hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

  MP_list_ref <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv")
  MP_list_ref <- as.list(MP_list_ref)
  mp_term2gene <- data.frame(
    term = rep(names(MP_list_ref), lengths(MP_list_ref)),
    gene = unlist(MP_list_ref),
    row.names = NULL
  )
  mp_term2gene$term <- sub("^MP", "3CA_mp", mp_term2gene$term)
  mp_term2name <- data.frame(
    term = unique(mp_term2gene$term),
    name = unique(mp_term2gene$term)
  )

  individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
  custom_files <- list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE)
  custom_refs <- lapply(custom_files, readRDS)
  names(custom_refs) <- sub(".*enrich_dev_", "", basename(custom_files)) %>% sub("\\.rds$", "", .)

  cluster_enrich <- lapply(names(mp_list), function(mp_name) {
    genes <- mp_list[[mp_name]]
    message(paste0("Processing MP: ", mp_name))

    res_GO <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                       ont = "BP", qvalueCutoff = 0.05, readable = TRUE)
    res_H <- enricher(gene = genes, TERM2GENE = hallmark_term2gene,
                      TERM2NAME = hallmark_term2name, qvalueCutoff = 0.05)
    res_M <- enricher(gene = genes, TERM2GENE = mp_term2gene,
                      TERM2NAME = mp_term2name, qvalueCutoff = 0.05)

    res_custom_list <- lapply(names(custom_refs), function(ref_name) {
      enricher(gene = genes, TERM2GENE = custom_refs[[ref_name]]$TERM2GENE,
               TERM2NAME = custom_refs[[ref_name]]$TERM2NAME,
               pAdjustMethod = "BH", qvalueCutoff = 0.05)
    })
    names(res_custom_list) <- names(custom_refs)

    c(list(rep_prog = mp_name, genes = genes, GO = res_GO, Hallmark = res_H, MPs_3CA = res_M),
      res_custom_list)
  })
  names(cluster_enrich) <- names(mp_list)

  enrich_heatmap <- function(cluster_enrich, element, top_per_program = 8, top_n = 80,
                              cap = 7, cols = viridis::magma(100, direction = -1),
                              fontsize_row = 7, fontsize_col = 9) {
    is_custom <- !element %in% c("GO", "Hallmark", "MPs_3CA")
    df_list <- lapply(names(cluster_enrich), function(prog) {
      er <- cluster_enrich[[prog]][[element]]
      if (is.null(er)) return(NULL)
      r <- tryCatch(er@result, error = function(e) NULL)
      if (is.null(r) || nrow(r) == 0) return(NULL)
      r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
      data_source <- if (is_custom) r else r_sig
      if (nrow(data_source) == 0 && !is_custom) return(NULL)
      term <- if ("Description" %in% colnames(data_source)) data_source$Description else data_source$ID
      data.frame(Program = prog, Term = term, padj = data_source$p.adjust,
                 Overlap = data_source$GeneRatio, stringsAsFactors = FALSE)
    })

    df <- dplyr::bind_rows(df_list)
    if (is.null(df) || nrow(df) == 0) {
      df <- data.frame(Program = character(), Term = character(),
                       padj = numeric(), Overlap = character(), stringsAsFactors = FALSE)
    }

    if (is_custom) {
      if (!element %in% names(custom_refs)) return(invisible(NULL))
      terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
    } else {
      if (nrow(df) == 0) return(invisible(NULL))
      terms_use <- df %>% dplyr::filter(padj < 0.05) %>%
        dplyr::arrange(Program, padj) %>% dplyr::group_by(Program) %>%
        dplyr::slice_head(n = top_per_program) %>% dplyr::ungroup() %>%
        dplyr::distinct(Term) %>% dplyr::pull(Term)
      if (length(terms_use) > top_n) {
        terms_use <- df %>% dplyr::filter(Term %in% terms_use) %>%
          dplyr::group_by(Term) %>%
          dplyr::summarise(min_p = min(padj), .groups = "drop") %>%
          dplyr::arrange(min_p) %>% dplyr::slice_head(n = top_n) %>%
          dplyr::pull(Term)
      }
    }

    if (!is.null(mp_tree_order)) {
      ordered_mps <- paste0("MP", mp_tree_order)
    } else {
      ordered_mps <- names(mp_list)
    }
    ordered_mps <- ordered_mps[ordered_mps %in% names(cluster_enrich)]

    full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)
    final_df <- full_grid %>%
      dplyr::left_join(df, by = c("Term", "Program")) %>%
      dplyr::mutate(
        score = tidyr::replace_na(pmin(-log10(padj), cap), 0),
        display_text = if (element %in% c("Hallmark", "GO", "MPs_3CA") || is_custom)
          tidyr::replace_na(Overlap, "") else ""
      )

    mat <- final_df %>% dplyr::select(Term, Program, score) %>%
      tidyr::pivot_wider(names_from = Program, values_from = score) %>%
      as.data.frame() %>% { row.names(.) <- .$Term; . } %>%
      dplyr::select(-Term) %>% as.matrix()

    text_mat <- final_df %>% dplyr::select(Term, Program, display_text) %>%
      tidyr::pivot_wider(names_from = Program, values_from = display_text) %>%
      as.data.frame() %>% { row.names(.) <- .$Term; . } %>%
      dplyr::select(-Term) %>% as.matrix()

    mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
    text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]
    if (nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))
    mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

    mp_sizes <- sapply(colnames(mat), function(x) length(mp_list[[x]]))
    col_labels <- paste0(colnames(mat), "\nn=", mp_sizes)

    cluster_rows_param <- FALSE; row_gaps <- NULL
    if (is_custom) {
      mat <- mat[terms_use, , drop = FALSE]
      text_mat <- text_mat[terms_use, , drop = FALSE]
    } else {
      best_mp <- colnames(mat)[max.col(mat, ties.method = "first")]
      row_order <- order(match(best_mp, colnames(mat)), -rowSums(mat))
      mat <- mat[row_order, , drop = FALSE]
      text_mat <- text_mat[row_order, , drop = FALSE]
      groups <- colnames(mat)[max.col(mat, ties.method = "first")]
      row_gaps <- which(groups[-length(groups)] != groups[-1])
    }

    breaks <- seq(0, cap, length.out = length(cols) + 1)
    pheatmap::pheatmap(mat, display_numbers = text_mat, number_color = "black",
                       fontsize_number = fontsize_row * 1.1, labels_col = col_labels,
                       color = cols, breaks = breaks, cluster_rows = cluster_rows_param,
                       cluster_cols = FALSE, gaps_row = row_gaps, border_color = NA,
                       show_colnames = TRUE, angle_col = 0,
                       fontsize_row = fontsize_row, fontsize_col = fontsize_col,
                       main = paste0(element, " Enrichment (-log10 padj)"))
    return(invisible(mat))
  }

  pdf(out_pdf, width = 12, height = 10)
  enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "GO",       top_per_program = 6, top_n = 60, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "MPs_3CA",  top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "Early_Embryogenesis", top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "Normal_Development_long", top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "Normal_Development_short", top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "Organogenesis_major", top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "Organogenesis_sub", top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "Adult_Epithelium", top_per_program = 8, top_n = 80, cols = cols_palette)
  enrich_heatmap(cluster_enrich, "Barretts_Oesophagus", top_per_program = 8, top_n = 80, cols = cols_palette)
  dev.off()
  cat("Saved combined PDF:", out_pdf, "\n")
}

cols_palette <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

####################
# Customized optimal-nMP similarity heatmap adapted from the current scRef
# centred workflow. PDO NMF programs are annotated by acquisition batch:
# treated/untreated samples are batch2, the Cynthia cohort is pdo, and the
# four new Souporcell samples are new4samples.
cat("Generating customized batch heatmap for optimal nMP:", optimal_nMP, "\n")
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(viridis)
})

sim_matrix <- geneNMF.metaprograms$programs.similarity
mp_clusters <- geneNMF.metaprograms$programs.clusters
keep_names <- names(mp_clusters)[!is.na(mp_clusters)]
ordered_names <- geneNMF.metaprograms$programs.tree$labels[
  geneNMF.metaprograms$programs.tree$order
]
final_ordered_names <- ordered_names[ordered_names %in% keep_names]
if (length(final_ordered_names) == 0) {
  stop("No NMF programs remained for the customized optimal-nMP heatmap.")
}
sim_matrix <- sim_matrix[final_ordered_names, final_ordered_names, drop = FALSE]

program_sample <- sub("\\.k[0-9]+\\.[0-9]+$", "", final_ordered_names)
program_batch <- ifelse(
  grepl("^TEMP_new4samples_|^SUR(1346|1363|1384|1391)(_|$)", program_sample),
  "new4samples",
  ifelse(grepl("_(Treated|Untreated)_PDO$", program_sample), "batch2", "pdo")
)

annotation_df <- data.frame(
  Program = final_ordered_names,
  Sample = program_sample,
  Metaprogram = paste0("MP", mp_clusters[final_ordered_names]),
  batch = factor(program_batch, levels = c("batch2", "pdo", "new4samples")),
  row.names = final_ordered_names,
  stringsAsFactors = FALSE
)
if (anyNA(annotation_df$batch)) {
  stop("At least one NMF program could not be assigned to a PDO batch.")
}
annotation_df$Metaprogram <- factor(
  annotation_df$Metaprogram,
  levels = unique(annotation_df$Metaprogram)
)

annotation_csv <- file.path(
  table_dir,
  paste0("Auto_centred_nMP_", optimal_nMP, "_program_batch_annotations.csv")
)
write.csv(annotation_df, annotation_csv, row.names = FALSE)

mp_cols <- setNames(
  colorRampPalette(brewer.pal(8, "Paired"))(length(levels(annotation_df$Metaprogram))),
  levels(annotation_df$Metaprogram)
)
batch_levels <- levels(annotation_df$batch)
batch_cols <- setNames(
  viridis::viridis(length(batch_levels), option = "turbo"),
  batch_levels
)

top_ha <- HeatmapAnnotation(
  df = annotation_df[, c("Metaprogram", "batch"), drop = FALSE],
  col = list(Metaprogram = mp_cols, batch = batch_cols),
  show_annotation_name = FALSE,
  show_legend = TRUE,
  simple_anno_size = grid::unit(2, "mm")
)
left_ha <- rowAnnotation(
  df = annotation_df[, c("Metaprogram", "batch"), drop = FALSE],
  col = list(Metaprogram = mp_cols, batch = batch_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  simple_anno_size = grid::unit(2, "mm")
)
col_fun <- colorRamp2(
  c(0.00, 0.12, 0.22, 0.70, 1.00),
  c("#FFFFFF", "#F6E8A6", "#E76F51", "#5E2A84", "#000000")
)

ht <- Heatmap(
  sim_matrix,
  name = "Similarity",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = annotation_df$Metaprogram,
  column_split = annotation_df$Metaprogram,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  rect_gp = grid::gpar(col = NA),
  border = FALSE,
  row_gap = grid::unit(0.4, "mm"),
  column_gap = grid::unit(0.4, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = top_ha,
  left_annotation = left_ha,
  use_raster = TRUE,
  raster_quality = 3,
  width = grid::unit(16, "cm"),
  height = grid::unit(16, "cm"),
  column_title_rot = 90,
  row_title_rot = 0,
  row_title_gp = grid::gpar(fontsize = 11),
  column_title_gp = grid::gpar(fontsize = 11)
)

custom_heatmap_pdf <- file.path(
  fig_dir,
  paste0("Auto_centred_nMP_", optimal_nMP, "_custom_heatmap_by_batch.pdf")
)
pdf(custom_heatmap_pdf, width = 10, height = 10)
draw(ht)
dev.off()
cat("Saved customized batch heatmap:", custom_heatmap_pdf, "\n")
cat("Saved program batch annotations:", annotation_csv, "\n")
####################

# Filter by silhouette < 0
mp_gene_lists <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  bad_mp_names <- paste0("MP", bad_mps)
  mp_gene_lists <- mp_gene_lists[!names(mp_gene_lists) %in% bad_mp_names]
}
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_gene_lists)))

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]

out_pdf <- file.path(fig_dir, "centred_initial_enrichment_anno.pdf")
run_enrichment_and_plot(mp_gene_lists, valid_cluster_ids, mp_tree_order, out_pdf, cols_palette)

cat("Auto_02_nmf_rank_selection_diagnostics.R completed successfully.\n")
