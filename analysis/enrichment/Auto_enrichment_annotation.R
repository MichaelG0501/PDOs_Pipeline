####################
# Enrichment annotation for PDO metaprograms
# Adapted from scRef_Pipeline/analysis/enrichment/enrichment_annotation.R
# Uses: Hallmark, GO:BP, 3CA MPs, developmental stage references
# Requires: optimal nMP result from find_optimal_nmf.R
####################
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

# ========================================
# 1. Load reference gene sets
# ========================================

# Hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

# 3CA metaprograms (cancer cell atlas)
MP_list <- read.csv("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv")
MP_list <- as.list(MP_list)
mp_term2gene <- data.frame(
  term = rep(names(MP_list), lengths(MP_list)),
  gene = unlist(MP_list),
  row.names = NULL
)
mp_term2gene$term <- sub("^MP", "3CA_mp", mp_term2gene$term)
mp_term2name <- data.frame(
  term = unique(mp_term2gene$term),
  name = unique(mp_term2gene$term)
)

# Developmental stage references
individual_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/"
custom_files <- list.files(individual_dir, pattern = "\\.rds$", full.names = TRUE)
custom_refs <- lapply(custom_files, readRDS)
names(custom_refs) <- sub(".*enrich_dev_", "", basename(custom_files)) %>% sub("\\.rds$", "", .)

# ========================================
# 2. Load metaprogram results
# ========================================

# Uses optimal nMP=13 result (determined by inflection-point analysis)
geneNMF.metaprograms <- readRDS("MP_outs_default.rds")
mp_gene_lists <- geneNMF.metaprograms$metaprograms.genes
mp_assignments <- geneNMF.metaprograms$programs.clusters

# Filter out MPs with silhouette < 0
print(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
mp_gene_lists <- mp_gene_lists[!names(mp_gene_lists) %in% bad_mp_names]

# Also filter mp_assignments to only keep valid cluster IDs
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_gene_lists)))
mp_assignments <- mp_assignments[mp_assignments %in% valid_cluster_ids & !is.na(mp_assignments)]

####################
# Additional filter: remove low-coverage MPs (<25% samples, i.e. <5/20)
####################
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
mp_gene_lists <- mp_gene_lists[!names(mp_gene_lists) %in% low_coverage_mps]

####################
# Keep mp_assignments aligned after low-coverage filtering
####################
valid_cluster_ids <- as.numeric(gsub("\\D", "", names(mp_gene_lists)))
mp_assignments <- mp_assignments[mp_assignments %in% valid_cluster_ids & !is.na(mp_assignments)]

message(paste0("Using ", length(mp_gene_lists), " metaprograms after silhouette filtering"))
message(paste0("Filtered out: ", paste(bad_mp_names, collapse = ", ")))

####################
# Report low-coverage filtering
####################
message(paste0("Coverage-filtered out (<25% samples): ", paste(low_coverage_mps, collapse = ", ")))
message(paste0("Using ", length(mp_gene_lists), " metaprograms after silhouette + coverage filtering"))

####################
# PDO MP display labels provided by user
####################
# mp_descriptions <- c(
#   "MP6"  = "MP6_G2M_mitotic",
#   "MP7"  = "MP7_DNA",
#   "MP5"  = "MP5_MYC Biosynth",
#   "MP1"  = "MP1_G2M_checkpoint",
#   "MP3"  = "MP3_G1S_Cycle",
#   "MP8"  = "MP8_Columnar progenitor",
#   "MP10" = "MP10_Stress-induced plasticity",
#   "MP9"  = "MP9_EMT_related",
#   "MP4"  = "MP4_Intest diff"
# )
mp_descriptions <- names(mp_gene_lists)
names(mp_descriptions) <- mp_descriptions

# ========================================
# 3. Enrichment analysis per metaprogram
# ========================================

cluster_enrich <- lapply(names(mp_gene_lists), function(mp_name) {

  genes <- mp_gene_lists[[mp_name]]
  mp_id <- as.numeric(gsub("\\D", "", mp_name))
  members <- names(mp_assignments)[mp_assignments == mp_id]

  message(paste0("Processing MP: ", mp_name))

  # Standard Enrichments
  res_GO <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                     ont = "BP", qvalueCutoff = 0.05, readable = TRUE)

  res_H <- enricher(gene = genes, TERM2GENE = hallmark_term2gene,
                    TERM2NAME = hallmark_term2name, qvalueCutoff = 0.05)

  res_M <- enricher(gene = genes, TERM2GENE = mp_term2gene,
                    TERM2NAME = mp_term2name, qvalueCutoff = 0.05)

  # Custom Marker Enrichment Loop (developmental stages)
  res_custom_list <- lapply(names(custom_refs), function(ref_name) {
    message(paste0("  -> Running custom enrichment: ", ref_name))

    enricher(
      gene          = genes,
      TERM2GENE     = custom_refs[[ref_name]]$TERM2GENE,
      TERM2NAME     = custom_refs[[ref_name]]$TERM2NAME,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05
    )
  })
  names(res_custom_list) <- names(custom_refs)

  # Consolidate all results into one flat list
  base_results <- list(
    rep_prog = mp_name,
    members  = members,
    genes    = genes,
    GO       = res_GO,
    Hallmark = res_H,
    MPs_3CA  = res_M
  )

  return(c(base_results, res_custom_list))
})

names(cluster_enrich) <- names(mp_gene_lists)
saveRDS(cluster_enrich, "cluster_enrich.rds")

# ========================================
# 4. Heatmap visualization function
# ========================================

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]

####################
# PDO order requirement: reverse tree order
####################
mp_tree_order <- rev(mp_tree_order)

enrich_heatmap <- function(cluster_enrich, element,
                           top_per_program = 8, top_n = 80, cap = 7,
                           cols = viridis::magma(100, direction = -1),
                           fontsize_row = 7, fontsize_col = 9) {

  is_custom <- !element %in% c("GO", "Hallmark", "MPs_3CA")

  df_list <- lapply(names(cluster_enrich), function(prog) {
    er <- cluster_enrich[[prog]][[element]]
    if (is.null(er)) return(NULL)

    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)

    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    data_source <- if(is_custom) r else r_sig

    if (nrow(data_source) == 0 && !is_custom) return(NULL)

    term <- if ("Description" %in% colnames(data_source)) data_source$Description else data_source$ID
    data.frame(
      Program = prog,
      Term = term,
      padj = data_source$p.adjust,
      Overlap = data_source$GeneRatio,
      stringsAsFactors = FALSE
    )
  })

  df <- dplyr::bind_rows(df_list)

  if (is.null(df) || nrow(df) == 0) {
    df <- data.frame(
      Program = character(),
      Term = character(),
      padj = numeric(),
      Overlap = character(),
      stringsAsFactors = FALSE
    )
  }

  # Define terms_use
  if (is_custom) {
    if (!element %in% names(custom_refs)) {
      message("Custom reference not found for element: ", element)
      return(invisible(NULL))
    }
    terms_use <- as.character(custom_refs[[element]]$TERM2NAME$term)
  } else {
    if (nrow(df) == 0) {
      message("No significant results found for: ", element)
      return(invisible(NULL))
    }
    terms_use <- df %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::arrange(Program, padj) %>%
      dplyr::group_by(Program) %>%
      dplyr::slice_head(n = top_per_program) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(Term) %>%
      dplyr::pull(Term)

    if (length(terms_use) > top_n) {
      terms_use <- df %>%
        dplyr::filter(Term %in% terms_use) %>%
        dplyr::group_by(Term) %>%
        dplyr::summarise(min_p = min(padj), .groups = "drop") %>%
        dplyr::arrange(min_p) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::pull(Term)
    }
  }

  # Matrix Construction
  ordered_mps <- paste0("MP", mp_tree_order)
  ordered_mps <- ordered_mps[ordered_mps %in% names(cluster_enrich)]
  full_grid <- expand.grid(Term = terms_use, Program = ordered_mps, stringsAsFactors = FALSE)

  final_df <- full_grid %>%
    dplyr::left_join(df, by = c("Term", "Program")) %>%
    dplyr::mutate(
      score = tidyr::replace_na(pmin(-log10(padj), cap), 0),
      display_text = if(element %in% c("Hallmark","GO","MPs_3CA") || is_custom) tidyr::replace_na(Overlap, "") else ""
    )

  # Score Matrix
  mat <- final_df %>%
    dplyr::select(Term, Program, score) %>%
    tidyr::pivot_wider(names_from = Program, values_from = score) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()

  # Display Text Matrix
  text_mat <- final_df %>%
    dplyr::select(Term, Program, display_text) %>%
    tidyr::pivot_wider(names_from = Program, values_from = display_text) %>%
    as.data.frame() %>% { row.names(.) <- .$Term; . } %>% dplyr::select(-Term) %>% as.matrix()

  # Align columns and rows
  mat <- mat[terms_use, ordered_mps[ordered_mps %in% colnames(mat)], drop = FALSE]
  text_mat <- text_mat[terms_use, colnames(mat), drop = FALSE]

  ####################
  # Guard: ensure numeric matrix for pheatmap scaling
  ####################
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message("No matrix content for element: ", element)
    return(invisible(NULL))
  }

  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))

  mp_sizes <- sapply(colnames(mat), function(x) {
    length(mp_gene_lists[[x]])
  })
  ####################
  # Use descriptive MP labels where available
  ####################
  display_mp <- ifelse(colnames(mat) %in% names(mp_descriptions), mp_descriptions[colnames(mat)], colnames(mat))
  col_labels <- paste0(display_mp, "\nn=", mp_sizes)

  # Sorting and Gaps
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

  # Render
  breaks <- seq(0, cap, length.out = length(cols) + 1)
  pheatmap::pheatmap(mat,
                     display_numbers = text_mat,
                     number_color = "black",
                     fontsize_number = fontsize_row * 1.1,
                     labels_col = col_labels,
                     color = cols,
                     breaks = breaks,
                     cluster_rows = cluster_rows_param,
                     cluster_cols = FALSE,
                     gaps_row = row_gaps,
                     border_color = NA,
                     show_colnames = TRUE,
                     angle_col = 0,
                     fontsize_row = fontsize_row,
                     fontsize_col = fontsize_col,
                     main = paste0(element, " Enrichment (-log10 padj)"))

  return(invisible(mat))
}

# ========================================
# 5. Generate enrichment heatmaps
# ========================================

cols_palette <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)

# Combined PDF
pdf("PDOs_enrichment_annotation.pdf", width = 10, height = 8)
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
cat("Saved combined PDF: PDOs_enrichment_annotation.pdf\n")

# Individual PNGs
png("enrich_Hallmark.png", width = 2000, height = 1750, res = 300)
enrich_heatmap(cluster_enrich, "Hallmark", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_GO.png", width = 2300, height = 2000, res = 300)
enrich_heatmap(cluster_enrich, "GO", top_per_program = 6, top_n = 60, cols = cols_palette)
dev.off()

png("enrich_MPs_3CA.png", width = 2000, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "MPs_3CA", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Early_Embryogenesis.png", width = 2200, height = 1500, res = 300)
enrich_heatmap(cluster_enrich, "Early_Embryogenesis", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Normal_Development_long.png", width = 2900, height = 3000, res = 300)
enrich_heatmap(cluster_enrich, "Normal_Development_long", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Normal_Development_short.png", width = 2900, height = 3000, res = 300)
enrich_heatmap(cluster_enrich, "Normal_Development_short", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Organogenesis_major.png", width = 2200, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "Organogenesis_major", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Organogenesis_sub.png", width = 2500, height = 1800, res = 300)
enrich_heatmap(cluster_enrich, "Organogenesis_sub", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Adult_Epithelium.png", width = 2250, height = 1900, res = 300)
enrich_heatmap(cluster_enrich, "Adult_Epithelium", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

png("enrich_Barretts_Oesophagus.png", width = 2250, height = 1900, res = 300)
enrich_heatmap(cluster_enrich, "Barretts_Oesophagus", top_per_program = 8, top_n = 80, cols = cols_palette)
dev.off()

message("All enrichment heatmaps saved to PDOs_outs/")
