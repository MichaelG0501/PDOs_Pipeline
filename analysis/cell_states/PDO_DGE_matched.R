####################
# Auto_PDO_DGE_matched.R
#
# Comprehensive DGE analysis for 4 matched Treated vs Untreated PDO sample pairs
#
# Three types of DGE analysis:
#   Type 1: Per patient - Treated vs Untreated (all cells)
#   Type 2: Per patient per state - Treated vs Untreated within each of 5 main states
#   Type 3: Per state - All patients combined - Treated vs Untreated (no overlap analysis)
#
# Outputs:
#   - Individual DEG lists per patient/state
#   - Overlap visualization (UpSet + Venn)
#   - Enrichment analysis (Hallmark + GO) for each list and overlaps
#   - Volcano plots and combined heatmaps
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(enrichplot)
library(pheatmap)
library(UpSetR)
library(VennDiagram)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(scales)
library(stringr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(viridis)
library(futile.logger)

# Suppress VennDiagram logging
flog.threshold(ERROR)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

message("=== Loading data ===")

# Load main Seurat object
pdos <- readRDS("PDOs_merged.rds")
state_vec <- readRDS("Auto_PDO_final_states.rds")
pdos$state <- state_vec[Cells(pdos)]

# Fix Batch
pdos$Batch <- ifelse(pdos$Batch %in% c("Treated_PDO", "Untreated_PDO"), "New_batch", "Cynthia_batch")

####################
# Configuration
####################

# Matched samples (4 pairs) - ordered by clinical response
matched_samples <- c(
  "SUR1070_Treated_PDO", "SUR1070_Untreated_PDO",
  "SUR1090_Treated_PDO", "SUR1090_Untreated_PDO",
  "SUR1072_Treated_PDO", "SUR1072_Untreated_PDO",
  "SUR1181_Treated_PDO", "SUR1181_Untreated_PDO"
)

# Ordered by clinical response with labels
patients <- c("SUR1070", "SUR1090", "SUR1072", "SUR1181")
patient_labels <- c(
  SUR1070 = "SUR1070 (Pre_Responder)",
  SUR1090 = "SUR1090 (Post_Responder)",
  SUR1072 = "SUR1072 (Post_Nonresponder)",
  SUR1181 = "SUR1181 (Pre_Nonresponder)"
)

# 5 main states (excluding Unresolved and Hybrid)
main_states <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "3CA_EMT_and_Protein_maturation"
)

# State colors (matching PDO_finalize_states.R)
state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

# Patient colors (in clinical response order) - now keyed by full label
patient_cols <- c(
  "SUR1070 (Pre_Responder)" = "#4C78A8",
  "SUR1090 (Post_Responder)" = "#B07AA1",
  "SUR1072 (Post_Nonresponder)" = "#59A14F",
  "SUR1181 (Pre_Nonresponder)" = "#F28E2B"
)

# DGE thresholds
logfc_threshold <- 0.25
padj_threshold <- 0.05
min_pct <- 0.1

####################
# Setup output directory
####################

out_dir <- "DGE_matched_analysis"
if (!dir.exists(out_dir)) dir.create(out_dir)

####################
# Subset to matched samples
####################

message("=== Subsetting to matched samples ===")
pdos_matched <- subset(pdos, subset = orig.ident %in% matched_samples)
pdos_matched$Treatment <- ifelse(grepl("Treated", pdos_matched$orig.ident), "Treated", "Untreated")
pdos_matched$Patient <- str_extract(pdos_matched$orig.ident, "^SUR\\d+")

message(sprintf("Total matched cells: %d", ncol(pdos_matched)))
print(table(pdos_matched$Patient, pdos_matched$Treatment))

####################
# Enrichment Setup (following enrichment_annotation.R pattern)
####################

message("=== Setting up enrichment databases ===")

# Hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term2gene <- hallmark_sets[, c("gs_name", "gene_symbol")]
hallmark_term2name <- hallmark_sets[, c("gs_name", "gs_name")]

# Function to run enrichment
run_enrichment <- function(genes, name_prefix) {
  results <- list()
  
  if (length(genes) < 5) {
    message(sprintf("  [%s] Too few genes (%d) for enrichment", name_prefix, length(genes)))
    return(results)
  }
  
  # GO Biological Process
  res_GO <- tryCatch({
    enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
  }, error = function(e) NULL)
  results$GO <- res_GO
  
  # Hallmark
  res_H <- tryCatch({
    enricher(
      gene = genes,
      TERM2GENE = hallmark_term2gene,
      TERM2NAME = hallmark_term2name,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05
    )
  }, error = function(e) NULL)
  results$Hallmark <- res_H
  
  return(results)
}

####################
# Helper: Volcano plot
####################

make_volcano <- function(deg_df, title, top_n = 20) {
  if (is.null(deg_df) || nrow(deg_df) == 0) return(NULL)
  
  deg_df <- deg_df %>%
    mutate(
      sig = case_when(
        p_val_adj < padj_threshold & avg_log2FC > logfc_threshold ~ "Up",
        p_val_adj < padj_threshold & avg_log2FC < -logfc_threshold ~ "Down",
        TRUE ~ "NS"
      ),
      neglog10p = -log10(p_val_adj + 1e-300)
    )
  
  # Top genes for labeling
  top_up <- deg_df %>% dplyr::filter(sig == "Up") %>% dplyr::slice_max(avg_log2FC, n = top_n/2)
  top_down <- deg_df %>% dplyr::filter(sig == "Down") %>% dplyr::slice_min(avg_log2FC, n = top_n/2)
  top_genes <- bind_rows(top_up, top_down)
  
  p <- ggplot(deg_df, aes(avg_log2FC, neglog10p)) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Up" = "firebrick3", "Down" = "steelblue", "NS" = "grey70")) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 2.5,
      max.overlaps = 30,
      segment.size = 0.2
    ) +
    labs(
      title = title,
      x = "log2 Fold Change (Treated vs Untreated)",
      y = "-log10(adjusted p-value)",
      color = "Regulation"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10)
    )
  
  return(p)
}

####################
# Helper: UpSet plot from gene lists
####################

make_upset <- function(gene_lists, title) {
  # Filter out empty lists
  gene_lists <- gene_lists[sapply(gene_lists, length) > 0]
  
  if (length(gene_lists) < 2) {
    message(sprintf("  [%s] Need at least 2 non-empty gene sets for UpSet", title))
    return(NULL)
  }
  
  # Create binary matrix
  all_genes <- unique(unlist(gene_lists))
  if (length(all_genes) == 0) return(NULL)
  
  mat <- sapply(gene_lists, function(g) as.integer(all_genes %in% g))
  rownames(mat) <- all_genes
  
  up <- upset(
    fromList(gene_lists),
    nsets = length(gene_lists),
    order.by = "freq",
    decreasing = TRUE,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "DEG Count",
    text.scale = 1.3,
    main.bar.color = "steelblue",
    sets.bar.color = unname(patient_cols[names(gene_lists)])
  )
  print(up)
  grid.text(title, x = 0.5, y = 0.95, gp = gpar(fontsize = 14, fontface = "bold"))
  return(invisible(up))
}

####################
# Helper: 4-way Venn diagram
####################

make_venn_4way <- function(gene_lists, title, filename) {
  # Ensure 4 lists
  gene_lists <- gene_lists[1:min(4, length(gene_lists))]
  
  if (length(gene_lists) < 2) {
    message(sprintf("  [%s] Need at least 2 sets for Venn", title))
    return(NULL)
  }
  
  # Pad with empty lists if needed
  while (length(gene_lists) < 4) {
    gene_lists[[paste0("Empty", length(gene_lists) + 1)]] <- character(0)
  }
  
  cols <- patient_cols[names(gene_lists)]
  cols[is.na(cols)] <- "grey70"
  
  venn.diagram(
    x = gene_lists,
    category.names = names(gene_lists),
    filename = filename,
    output = TRUE,
    imagetype = "png",
    height = 2000,
    width = 2000,
    resolution = 300,
    lwd = 2,
    fill = unname(cols),
    alpha = 0.5,
    cex = 1.2,
    fontface = "bold",
    cat.cex = 1.0,
    cat.fontface = "bold",
    main = title,
    main.cex = 1.5
  )
}

####################
# Helper: Enrichment heatmap (following enrichment_annotation.R)
####################

plot_enrich_heatmap <- function(enrich_list, element, title,
                                 top_per_set = 5, top_n = 40, cap = 7) {
  
  df_list <- lapply(names(enrich_list), function(set_name) {
    er <- enrich_list[[set_name]][[element]]
    if (is.null(er)) return(NULL)
    
    r <- tryCatch(er@result, error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    
    r_sig <- r[which(r$p.adjust < 0.05 & r$p.adjust > 0), ]
    if (nrow(r_sig) == 0) return(NULL)
    
    term <- if ("Description" %in% colnames(r_sig)) r_sig$Description else r_sig$ID
    data.frame(
      Set = set_name,
      Term = term,
      padj = r_sig$p.adjust,
      GeneRatio = r_sig$GeneRatio,
      stringsAsFactors = FALSE
    )
  })
  
  df <- bind_rows(df_list)
  if (nrow(df) == 0) {
    message(sprintf("  No significant %s results for heatmap", element))
    return(NULL)
  }
  
  # Select top terms
  terms_use <- df %>%
    dplyr::filter(padj < 0.05) %>%
    arrange(Set, padj) %>%
    group_by(Set) %>%
    dplyr::slice_head(n = top_per_set) %>%
    ungroup() %>%
    distinct(Term) %>%
    pull(Term)
  
  if (length(terms_use) > top_n) {
    terms_use <- df %>%
      dplyr::filter(Term %in% terms_use) %>%
      group_by(Term) %>%
      summarise(min_p = min(padj), .groups = "drop") %>%
      arrange(min_p) %>%
      dplyr::slice_head(n = top_n) %>%
      pull(Term)
  }
  
  # Build matrix
  full_grid <- expand.grid(Term = terms_use, Set = names(enrich_list), stringsAsFactors = FALSE)
  final_df <- full_grid %>%
    left_join(df, by = c("Term", "Set")) %>%
    mutate(score = replace_na(pmin(-log10(padj), cap), 0))
  
  mat <- final_df %>%
    dplyr::select(Term, Set, score) %>%
    pivot_wider(names_from = Set, values_from = score) %>%
    column_to_rownames("Term") %>%
    as.matrix()
  
  if (nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
  
  # Sort rows by best column
  best_col <- colnames(mat)[max.col(mat, ties.method = "first")]
  row_order <- order(match(best_col, colnames(mat)), -rowSums(mat))
  mat <- mat[row_order, , drop = FALSE]
  
  cols_palette <- colorRampPalette(c("#ffffff", "#ffcccc", "#ff6666", "#cc0000", "#660000"))(100)
  
  pheatmap::pheatmap(
    mat,
    color = cols_palette,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = 10,
    angle_col = 45,
    main = paste0(title, " - ", element, " Enrichment")
  )
}

####################
# TYPE 1: DGE - Treated vs Untreated per patient (all cells)
####################

message("\n=== TYPE 1: DGE - Treated vs Untreated per patient ===")

dge_type1_results <- list()
dge_type1_up <- list()
dge_type1_down <- list()

for (pat in patients) {
  message(sprintf("\nProcessing patient: %s", pat))
  
  # Subset to this patient
  pat_cells <- Cells(pdos_matched)[pdos_matched$Patient == pat]
  pat_obj <- subset(pdos_matched, cells = pat_cells)
  
  n_treated <- sum(pat_obj$Treatment == "Treated")
  n_untreated <- sum(pat_obj$Treatment == "Untreated")
  message(sprintf("  Treated: %d cells, Untreated: %d cells", n_treated, n_untreated))
  
  if (n_treated < 30 || n_untreated < 30) {
    message(sprintf("  [%s] Insufficient cells, skipping", pat))
    next
  }
  
  # Set identity
  Idents(pat_obj) <- "Treatment"
  DefaultAssay(pat_obj) <- "RNA"
  
  # Find variable features
  pat_obj <- FindVariableFeatures(pat_obj, nfeatures = 3000, verbose = FALSE)
  
  # Run DGE: Treated vs Untreated
  dge <- tryCatch({
    FindMarkers(
      pat_obj,
      ident.1 = "Treated",
      ident.2 = "Untreated",
      test.use = "wilcox",
      logfc.threshold = 0,  # Get all for volcano
      min.pct = min_pct,
      features = VariableFeatures(pat_obj),
      verbose = FALSE
    )
  }, error = function(e) {
    message(sprintf("  [%s] FindMarkers error: %s", pat, e$message))
    NULL
  })
  
  if (is.null(dge) || nrow(dge) == 0) {
    message(sprintf("  [%s] No DEGs found", pat))
    next
  }
  
  pat_label <- patient_labels[pat]
  dge$gene <- rownames(dge)
  dge$patient <- pat_label
  dge_type1_results[[pat_label]] <- dge
  
  # Significant DEGs (upregulated means upregulated in Treated)
  sig_up <- dge %>% dplyr::filter(p_val_adj < padj_threshold, avg_log2FC > logfc_threshold)
  sig_down <- dge %>% dplyr::filter(p_val_adj < padj_threshold, avg_log2FC < -logfc_threshold)
  
  dge_type1_up[[pat_label]] <- sig_up$gene
  dge_type1_down[[pat_label]] <- sig_down$gene
  
  message(sprintf("  [%s] Up: %d, Down: %d DEGs", pat_label, nrow(sig_up), nrow(sig_down)))
}

# Save Type 1 results
saveRDS(dge_type1_results, file.path(out_dir, "DGE_Type1_results.rds"))

####################
# TYPE 1: Volcano plots
####################

message("\n=== TYPE 1: Generating volcano plots ===")

pdf(file.path(out_dir, "Type1_volcano_plots.pdf"), width = 12, height = 10)

volcano_plots <- list()
for (pat in names(dge_type1_results)) {
  p <- make_volcano(dge_type1_results[[pat]], paste0(pat, ": Treated vs Untreated"))
  if (!is.null(p)) volcano_plots[[pat]] <- p
}

if (length(volcano_plots) > 0) {
  combined <- wrap_plots(volcano_plots, ncol = 2) + 
    plot_annotation(title = "Type 1 DGE: Treated vs Untreated per Patient",
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14)))
  print(combined)
}
dev.off()

####################
# TYPE 1: Overlap visualization
####################

message("\n=== TYPE 1: Overlap visualization ===")

# UpSet for upregulated genes
pdf(file.path(out_dir, "Type1_UpSet_up.pdf"), width = 10, height = 6)
if (sum(sapply(dge_type1_up, length) > 0) >= 2) {
  make_upset(dge_type1_up, "Upregulated DEGs")
}
dev.off()

# UpSet for downregulated genes
pdf(file.path(out_dir, "Type1_UpSet_down.pdf"), width = 10, height = 6)
if (sum(sapply(dge_type1_down, length) > 0) >= 2) {
  make_upset(dge_type1_down, "Downregulated DEGs")
}
dev.off()

# Venn diagrams
if (sum(sapply(dge_type1_up, length) > 0) >= 2) {
  make_venn_4way(dge_type1_up, "Upregulated DEGs (Treated vs Untreated)",
                 file.path(out_dir, "Type1_Venn_up.png"))
}

if (sum(sapply(dge_type1_down, length) > 0) >= 2) {
  make_venn_4way(dge_type1_down, "Downregulated DEGs (Treated vs Untreated)",
                 file.path(out_dir, "Type1_Venn_down.png"))
}

####################
# TYPE 1: Enrichment analysis
####################

message("\n=== TYPE 1: Enrichment analysis ===")

# Enrichment for each patient
enrich_type1_up <- list()
enrich_type1_down <- list()

for (pat_label in names(dge_type1_up)) {
  message(sprintf("  Enrichment for %s (Up)", pat_label))
  enrich_type1_up[[pat_label]] <- run_enrichment(dge_type1_up[[pat_label]], paste0(pat_label, "_up"))
  
  message(sprintf("  Enrichment for %s (Down)", pat_label))
  enrich_type1_down[[pat_label]] <- run_enrichment(dge_type1_down[[pat_label]], paste0(pat_label, "_down"))
}

# Find overlapping genes in 3+ patients (instead of all)
find_overlap_3plus <- function(gene_lists) {
  gene_lists <- gene_lists[sapply(gene_lists, length) > 0]
  if (length(gene_lists) < 3) return(character(0))
  all_genes <- unique(unlist(gene_lists))
  gene_counts <- sapply(all_genes, function(g) sum(sapply(gene_lists, function(l) g %in% l)))
  all_genes[gene_counts >= 3]
}

type1_overlap_up <- find_overlap_3plus(dge_type1_up)
type1_overlap_down <- find_overlap_3plus(dge_type1_down)

message(sprintf("  Overlapping upregulated genes (3+ patients): %d", length(type1_overlap_up)))
message(sprintf("  Overlapping downregulated genes (3+ patients): %d", length(type1_overlap_down)))

# Enrichment for overlaps
enrich_type1_up[["Overlap_3plus"]] <- run_enrichment(type1_overlap_up, "Overlap_3plus_up")
enrich_type1_down[["Overlap_3plus"]] <- run_enrichment(type1_overlap_down, "Overlap_3plus_down")

# Save enrichment results
saveRDS(list(up = enrich_type1_up, down = enrich_type1_down), 
        file.path(out_dir, "DGE_Type1_enrichment.rds"))

# Enrichment heatmaps
pdf(file.path(out_dir, "Type1_enrichment_heatmaps.pdf"), width = 14, height = 10)
plot_enrich_heatmap(enrich_type1_up, "Hallmark", "Type 1 Upregulated")
plot_enrich_heatmap(enrich_type1_up, "GO", "Type 1 Upregulated")
plot_enrich_heatmap(enrich_type1_down, "Hallmark", "Type 1 Downregulated")
plot_enrich_heatmap(enrich_type1_down, "GO", "Type 1 Downregulated")
dev.off()

####################
# TYPE 2: DGE - Treated vs Untreated within each State per patient
####################

message("\n=== TYPE 2: DGE - Treated vs Untreated within each State per patient ===")

dge_type2_results <- list()
dge_type2_up <- list()
dge_type2_down <- list()

for (state in main_states) {
  message(sprintf("\n--- State: %s ---", state))
  
  state_key <- gsub(" ", "_", gsub("\\.", "", state))
  dge_type2_results[[state_key]] <- list()
  dge_type2_up[[state_key]] <- list()
  dge_type2_down[[state_key]] <- list()
  
  for (pat in patients) {
    message(sprintf("  Patient: %s", pat))
    
    # Subset to this patient and state
    pat_state_cells <- Cells(pdos_matched)[
      pdos_matched$Patient == pat & pdos_matched$state == state
    ]
    
    if (length(pat_state_cells) < 60) {
      message(sprintf("    [%s|%s] Insufficient cells (%d), skipping", pat, state, length(pat_state_cells)))
      next
    }
    
    pat_obj <- subset(pdos_matched, cells = pat_state_cells)
    
    n_treated <- sum(pat_obj$Treatment == "Treated")
    n_untreated <- sum(pat_obj$Treatment == "Untreated")
    message(sprintf("    Treated: %d, Untreated: %d", n_treated, n_untreated))
    
    if (n_treated < 20 || n_untreated < 20) {
      message(sprintf("    [%s|%s] Insufficient per-treatment cells, skipping", pat, state))
      next
    }
    
    # Set identity
    Idents(pat_obj) <- "Treatment"
    DefaultAssay(pat_obj) <- "RNA"
    
    # Find variable features
    pat_obj <- FindVariableFeatures(pat_obj, nfeatures = 3000, verbose = FALSE)
    
    # Run DGE
    dge <- tryCatch({
      FindMarkers(
        pat_obj,
        ident.1 = "Treated",
        ident.2 = "Untreated",
        test.use = "wilcox",
        logfc.threshold = 0,
        min.pct = min_pct,
        features = VariableFeatures(pat_obj),
        verbose = FALSE
      )
    }, error = function(e) {
      message(sprintf("    [%s|%s] FindMarkers error: %s", pat, state, e$message))
      NULL
    })
    
    if (is.null(dge) || nrow(dge) == 0) {
      message(sprintf("    [%s|%s] No DEGs found", pat, state))
      next
    }
    
    pat_label <- patient_labels[pat]
    dge$gene <- rownames(dge)
    dge$patient <- pat_label
    dge$state <- state
    dge_type2_results[[state_key]][[pat_label]] <- dge
    
    # Significant DEGs (upregulated means upregulated in Treated)
    sig_up <- dge %>% dplyr::filter(p_val_adj < padj_threshold, avg_log2FC > logfc_threshold)
    sig_down <- dge %>% dplyr::filter(p_val_adj < padj_threshold, avg_log2FC < -logfc_threshold)
    
    dge_type2_up[[state_key]][[pat_label]] <- sig_up$gene
    dge_type2_down[[state_key]][[pat_label]] <- sig_down$gene
    
    message(sprintf("    [%s|%s] Up: %d, Down: %d DEGs", pat_label, state, nrow(sig_up), nrow(sig_down)))
  }
}

# Save Type 2 results
saveRDS(dge_type2_results, file.path(out_dir, "DGE_Type2_results.rds"))

####################
# TYPE 2: Volcano plots per state
####################

message("\n=== TYPE 2: Generating volcano plots per state ===")

pdf(file.path(out_dir, "Type2_volcano_plots_by_state.pdf"), width = 14, height = 12)

for (state_key in names(dge_type2_results)) {
  state_results <- dge_type2_results[[state_key]]
  if (length(state_results) == 0) next
  
  volcano_plots <- list()
  for (pat in names(state_results)) {
    p <- make_volcano(state_results[[pat]], 
                      paste0(pat, "\n", gsub("_", " ", state_key)))
    if (!is.null(p)) volcano_plots[[pat]] <- p
  }
  
  if (length(volcano_plots) > 0) {
    combined <- wrap_plots(volcano_plots, ncol = 2) +
      plot_annotation(
        title = paste0("Type 2 DGE: ", gsub("_", " ", state_key)),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
      )
    print(combined)
  }
}
dev.off()

####################
# TYPE 2: Overlap visualization per state
####################

message("\n=== TYPE 2: Overlap visualization per state ===")

pdf(file.path(out_dir, "Type2_UpSet_by_state.pdf"), width = 10, height = 6)
for (state_key in names(dge_type2_up)) {
  up_lists <- dge_type2_up[[state_key]]
  down_lists <- dge_type2_down[[state_key]]
  
  if (sum(sapply(up_lists, length) > 0) >= 2) {
    message(sprintf("  UpSet for %s (Up)", state_key))
    make_upset(up_lists, paste0(gsub("_", " ", state_key), " - Upregulated"))
  }
  
  if (sum(sapply(down_lists, length) > 0) >= 2) {
    message(sprintf("  UpSet for %s (Down)", state_key))
    make_upset(down_lists, paste0(gsub("_", " ", state_key), " - Downregulated"))
  }
}
dev.off()

# Venn diagrams per state
for (state_key in names(dge_type2_up)) {
  up_lists <- dge_type2_up[[state_key]]
  down_lists <- dge_type2_down[[state_key]]
  
  if (sum(sapply(up_lists, length) > 0) >= 2) {
    make_venn_4way(up_lists, paste0(gsub("_", " ", state_key), " - Upregulated"),
                   file.path(out_dir, paste0("Type2_Venn_up_", state_key, ".png")))
  }
  
  if (sum(sapply(down_lists, length) > 0) >= 2) {
    make_venn_4way(down_lists, paste0(gsub("_", " ", state_key), " - Downregulated"),
                   file.path(out_dir, paste0("Type2_Venn_down_", state_key, ".png")))
  }
}

####################
# TYPE 2: Enrichment analysis
####################

message("\n=== TYPE 2: Enrichment analysis per state ===")

enrich_type2_up <- list()
enrich_type2_down <- list()

for (state_key in names(dge_type2_up)) {
  message(sprintf("\n  --- State: %s ---", state_key))
  
  enrich_type2_up[[state_key]] <- list()
  enrich_type2_down[[state_key]] <- list()
  
  up_lists <- dge_type2_up[[state_key]]
  down_lists <- dge_type2_down[[state_key]]
  
  # Per patient enrichment (with clinical response labels)
  for (pat_label in names(up_lists)) {
    message(sprintf("    Enrichment for %s (Up)", pat_label))
    enrich_type2_up[[state_key]][[pat_label]] <- run_enrichment(up_lists[[pat_label]], paste0(state_key, "_", pat_label, "_up"))
    
    message(sprintf("    Enrichment for %s (Down)", pat_label))
    enrich_type2_down[[state_key]][[pat_label]] <- run_enrichment(down_lists[[pat_label]], paste0(state_key, "_", pat_label, "_down"))
  }
  
  # Overlap enrichment (3+ patients instead of all)
  if (sum(sapply(up_lists, length) > 0) >= 3) {
    overlap_up <- find_overlap_3plus(up_lists)
    message(sprintf("    Overlapping upregulated genes (3+ patients): %d", length(overlap_up)))
    enrich_type2_up[[state_key]][["Overlap_3plus"]] <- run_enrichment(overlap_up, paste0(state_key, "_Overlap_3plus_up"))
  }
  
  if (sum(sapply(down_lists, length) > 0) >= 3) {
    overlap_down <- find_overlap_3plus(down_lists)
    message(sprintf("    Overlapping downregulated genes (3+ patients): %d", length(overlap_down)))
    enrich_type2_down[[state_key]][["Overlap_3plus"]] <- run_enrichment(overlap_down, paste0(state_key, "_Overlap_3plus_down"))
  }
}

# Save enrichment results
saveRDS(list(up = enrich_type2_up, down = enrich_type2_down),
        file.path(out_dir, "DGE_Type2_enrichment.rds"))

# Enrichment heatmaps per state
pdf(file.path(out_dir, "Type2_enrichment_heatmaps.pdf"), width = 14, height = 10)
for (state_key in names(enrich_type2_up)) {
  up_enrich <- enrich_type2_up[[state_key]]
  down_enrich <- enrich_type2_down[[state_key]]
  
  if (length(up_enrich) > 0) {
    plot_enrich_heatmap(up_enrich, "Hallmark", paste0("Type 2 ", gsub("_", " ", state_key), " - Upregulated"))
    plot_enrich_heatmap(up_enrich, "GO", paste0("Type 2 ", gsub("_", " ", state_key), " - Upregulated"))
  }
  
  if (length(down_enrich) > 0) {
    plot_enrich_heatmap(down_enrich, "Hallmark", paste0("Type 2 ", gsub("_", " ", state_key), " - Downregulated"))
    plot_enrich_heatmap(down_enrich, "GO", paste0("Type 2 ", gsub("_", " ", state_key), " - Downregulated"))
  }
}
dev.off()

####################
# Combined Heatmap: Top DEGs across patients
####################

message("\n=== Generating combined DEG heatmaps ===")

# Function to create combined expression heatmap
make_deg_heatmap <- function(deg_list, seurat_obj, title, top_n = 30) {
  # Get union of top DEGs
  all_degs <- lapply(deg_list, function(df) {
    if (is.null(df)) return(character(0))
    df %>%
      dplyr::filter(p_val_adj < padj_threshold, abs(avg_log2FC) > logfc_threshold) %>%
      dplyr::slice_max(abs(avg_log2FC), n = top_n) %>%
      pull(gene)
  })
  
  top_genes <- unique(unlist(all_degs))
  if (length(top_genes) == 0) return(NULL)
  
  # Get expression matrix
  expr_mat <- GetAssayData(seurat_obj, layer = "data")[top_genes, , drop = FALSE]
  
  # Average by patient x treatment
  meta <- seurat_obj@meta.data %>%
    mutate(group = paste(Patient, Treatment, sep = "_"))
  
  avg_expr <- sapply(unique(meta$group), function(g) {
    cells <- rownames(meta)[meta$group == g]
    rowMeans(as.matrix(expr_mat[, cells, drop = FALSE]))
  })
  
  # Z-score normalize
  avg_expr_z <- t(scale(t(avg_expr)))
  
  # Order columns: all Untreated then all Treated, following patient order
  col_order <- c(paste0(patients, "_Untreated"), paste0(patients, "_Treated"))
  col_order <- col_order[col_order %in% colnames(avg_expr_z)]
  avg_expr_z <- avg_expr_z[, col_order, drop = FALSE]
  
  # Update column names to include full patient labels
  new_colnames <- sapply(colnames(avg_expr_z), function(x) {
    p <- str_extract(x, "^SUR\\d+")
    t_val <- ifelse(grepl("Treated", x), "Treated", "Untreated")
    paste0(patient_labels[p], "_", t_val)
  })
  colnames(avg_expr_z) <- new_colnames
  
  # Column annotation (created AFTER reordering to match)
  col_df <- data.frame(
    Group = colnames(avg_expr_z),
    Patient = str_extract(colnames(avg_expr_z), "^SUR\\d+ \\([^\\)]+\\)"),
    Treatment = ifelse(grepl("Treated", colnames(avg_expr_z)), "Treated", "Untreated")
  )
  rownames(col_df) <- col_df$Group
  
  col_ann <- HeatmapAnnotation(
    Patient = col_df$Patient,
    Treatment = col_df$Treatment,
    col = list(
      Patient = patient_cols,
      Treatment = c(Treated = "#4D4D4D", Untreated = "#E69F00")
    ),
    show_annotation_name = TRUE
  )
  
  # Heatmap
  Heatmap(
    avg_expr_z,
    name = "Z-score",
    top_annotation = col_ann,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 45,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 9),
    col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
}

pdf(file.path(out_dir, "Combined_DEG_heatmaps.pdf"), width = 12, height = 14)

# Type 1 combined heatmap
if (length(dge_type1_results) > 0) {
  ht <- make_deg_heatmap(dge_type1_results, pdos_matched, 
                         "Type 1: Top DEGs across Patients (Treated vs Untreated)")
  if (!is.null(ht)) draw(ht)
}

# Type 2 combined heatmap per state
for (state_key in names(dge_type2_results)) {
  state_results <- dge_type2_results[[state_key]]
  if (length(state_results) == 0) next
  
  # Subset to state cells
  state_cells <- Cells(pdos_matched)[pdos_matched$state == gsub("_", " ", state_key) |
                                      pdos_matched$state == state_key]
  if (length(state_cells) < 100) next
  
  state_obj <- subset(pdos_matched, cells = state_cells)
  
  ht <- make_deg_heatmap(state_results, state_obj,
                         paste0("Type 2: Top DEGs - ", gsub("_", " ", state_key)))
  if (!is.null(ht)) draw(ht)
}

dev.off()

####################
# Summary statistics table
####################

message("\n=== Generating summary statistics ===")

# Type 1 summary
type1_summary <- data.frame(
  Patient = names(dge_type1_results),
  Total_DEGs = sapply(dge_type1_results, function(df) {
    sum(df$p_val_adj < padj_threshold & abs(df$avg_log2FC) > logfc_threshold)
  }),
  Upregulated = sapply(dge_type1_up, length),
  Downregulated = sapply(dge_type1_down, length)
)

# Enforce patient order in Type 1 summary
type1_summary$Patient <- factor(type1_summary$Patient, levels = patient_labels[patients])

# Type 2 summary
type2_summary_list <- list()
for (state_key in names(dge_type2_results)) {
  for (pat_label in names(dge_type2_results[[state_key]])) {
    df <- dge_type2_results[[state_key]][[pat_label]]
    type2_summary_list[[paste0(state_key, "_", pat_label)]] <- data.frame(
      State = gsub("_", " ", state_key),
      Patient = pat_label,
      Total_DEGs = sum(df$p_val_adj < padj_threshold & abs(df$avg_log2FC) > logfc_threshold),
      Upregulated = length(dge_type2_up[[state_key]][[pat_label]]),
      Downregulated = length(dge_type2_down[[state_key]][[pat_label]])
    )
  }
}
type2_summary <- bind_rows(type2_summary_list)

# Enforce patient order in Type 2 summary
if (nrow(type2_summary) > 0) {
  type2_summary$Patient <- factor(type2_summary$Patient, levels = patient_labels[patients])
}

# Save summaries
write.csv(type1_summary, file.path(out_dir, "Type1_summary.csv"), row.names = FALSE)
write.csv(type2_summary, file.path(out_dir, "Type2_summary.csv"), row.names = FALSE)

####################
# Summary visualization
####################

pdf(file.path(out_dir, "Summary_barplots.pdf"), width = 12, height = 8)

# Type 1 summary plot
p1 <- type1_summary %>%
  pivot_longer(cols = c(Upregulated, Downregulated), names_to = "Direction", values_to = "Count") %>%
  ggplot(aes(Patient, Count, fill = Direction)) +
    geom_col(position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c(Upregulated = "firebrick3", Downregulated = "steelblue")) +
  labs(title = "Type 1: DEG counts per patient (Treated vs Untreated)",
       x = NULL, y = "Number of DEGs") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
  )

print(p1)

# Type 2 summary plot
if (nrow(type2_summary) > 0) {
  p2 <- type2_summary %>%
    pivot_longer(cols = c(Upregulated, Downregulated), names_to = "Direction", values_to = "Count") %>%
    ggplot(aes(Patient, Count, fill = Direction)) +
    geom_col(position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c(Upregulated = "firebrick3", Downregulated = "steelblue")) +
    facet_wrap(~ State, scales = "free_y", ncol = 2) +
    labs(title = "Type 2: DEG counts per patient per state",
         x = NULL, y = "Number of DEGs") +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
    )
  
  print(p2)
}

dev.off()

####################
# TYPE 3: DGE - All patients combined per State (Treated vs Untreated)
# No overlap analysis - just expression heatmap and enrichment
####################

message("\n=== TYPE 3: DGE - All patients combined per State ===")

dge_type3_results <- list()
dge_type3_up <- list()
dge_type3_down <- list()

for (state in main_states) {
  message(sprintf("\n--- State: %s ---", state))
  
  state_key <- gsub(" ", "_", gsub("\\.", "", state))
  
  # Subset to this state across ALL patients
  state_cells <- Cells(pdos_matched)[pdos_matched$state == state]
  
  if (length(state_cells) < 100) {
    message(sprintf("  [%s] Insufficient cells (%d), skipping", state, length(state_cells)))
    next
  }
  
  state_obj <- subset(pdos_matched, cells = state_cells)
  
  n_treated <- sum(state_obj$Treatment == "Treated")
  n_untreated <- sum(state_obj$Treatment == "Untreated")
  message(sprintf("  Total: %d cells (Treated: %d, Untreated: %d)", 
                  length(state_cells), n_treated, n_untreated))
  
  if (n_treated < 50 || n_untreated < 50) {
    message(sprintf("  [%s] Insufficient per-treatment cells, skipping", state))
    next
  }
  
  # Set identity
  Idents(state_obj) <- "Treatment"
  DefaultAssay(state_obj) <- "RNA"
  
  # Find variable features
  state_obj <- FindVariableFeatures(state_obj, nfeatures = 3000, verbose = FALSE)
  
  # Run DGE: Treated vs Untreated
  dge <- tryCatch({
    FindMarkers(
      state_obj,
      ident.1 = "Treated",
      ident.2 = "Untreated",
      test.use = "wilcox",
      logfc.threshold = 0,  # Get all for heatmap
      min.pct = min_pct,
      features = VariableFeatures(state_obj),
      verbose = FALSE
    )
  }, error = function(e) {
    message(sprintf("  [%s] FindMarkers error: %s", state, e$message))
    NULL
  })
  
  if (is.null(dge) || nrow(dge) == 0) {
    message(sprintf("  [%s] No DEGs found", state))
    next
  }
  
  dge$gene <- rownames(dge)
  dge$state <- state
  dge_type3_results[[state_key]] <- dge
  
  # Significant DEGs
  sig_up <- dge %>% dplyr::filter(p_val_adj < padj_threshold, avg_log2FC > logfc_threshold)
  sig_down <- dge %>% dplyr::filter(p_val_adj < padj_threshold, avg_log2FC < -logfc_threshold)
  
  dge_type3_up[[state_key]] <- sig_up$gene
  dge_type3_down[[state_key]] <- sig_down$gene
  
  message(sprintf("  [%s] Up: %d, Down: %d DEGs", state, nrow(sig_up), nrow(sig_down)))
}

# Save Type 3 results
saveRDS(dge_type3_results, file.path(out_dir, "DGE_Type3_results.rds"))

####################
# TYPE 3: Expression Heatmap (top DEGs by significance, averaged by Treatment)
####################

message("\n=== TYPE 3: Generating expression heatmaps ===")

make_type3_heatmap <- function(dge_df, seurat_obj, state_name, top_n = 40) {
  if (is.null(dge_df) || nrow(dge_df) == 0) return(NULL)
  
  # Select top DEGs by significance (balanced up/down)
  sig_degs <- dge_df %>%
    dplyr::filter(p_val_adj < padj_threshold, abs(avg_log2FC) > logfc_threshold) %>%
    arrange(p_val_adj)
  
  if (nrow(sig_degs) == 0) return(NULL)
  
  # Take top N, balanced between up and down
  top_up <- sig_degs %>% dplyr::filter(avg_log2FC > 0) %>% dplyr::slice_head(n = top_n/2)
  top_down <- sig_degs %>% dplyr::filter(avg_log2FC < 0) %>% dplyr::slice_head(n = top_n/2)
  top_genes <- c(top_up$gene, top_down$gene)
  
  if (length(top_genes) == 0) return(NULL)
  
  # Get expression matrix
  expr_mat <- GetAssayData(seurat_obj, layer = "data")[top_genes, , drop = FALSE]
  
  # Average by Treatment only (2 columns)
  avg_expr <- sapply(c("Untreated", "Treated"), function(trt) {
    cells <- Cells(seurat_obj)[seurat_obj$Treatment == trt]
    rowMeans(as.matrix(expr_mat[, cells, drop = FALSE]))
  })
  
  # Z-score normalize
  avg_expr_z <- t(scale(t(avg_expr)))
  
  # Column annotation
  col_ann <- HeatmapAnnotation(
    Treatment = colnames(avg_expr_z),
    col = list(Treatment = c(Treated = "#4D4D4D", Untreated = "#E69F00")),
    show_annotation_name = TRUE
  )
  
  # Row annotation with direction
  row_direction <- ifelse(dge_df[top_genes, "avg_log2FC"] > 0, "Up", "Down")
  row_ann <- rowAnnotation(
    Direction = row_direction,
    col = list(Direction = c(Up = "firebrick3", Down = "steelblue")),
    show_annotation_name = TRUE
  )
  
  # Order rows by direction then by log2FC
  row_order <- c(
    top_genes[top_genes %in% top_up$gene][order(-dge_df[top_genes[top_genes %in% top_up$gene], "avg_log2FC"])],
    top_genes[top_genes %in% top_down$gene][order(dge_df[top_genes[top_genes %in% top_down$gene], "avg_log2FC"])]
  )
  avg_expr_z <- avg_expr_z[row_order, , drop = FALSE]
  row_direction <- row_direction[match(row_order, top_genes)]
  
  row_ann <- rowAnnotation(
    Direction = row_direction,
    col = list(Direction = c(Up = "firebrick3", Down = "steelblue")),
    show_annotation_name = TRUE
  )
  
  # Heatmap
  Heatmap(
    avg_expr_z,
    name = "Z-score",
    top_annotation = col_ann,
    right_annotation = row_ann,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 12),
    col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3")),
    column_title = paste0("Type 3: ", state_name, "\n(All patients combined)"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    width = unit(4, "cm")
  )
}

pdf(file.path(out_dir, "Type3_expression_heatmaps.pdf"), width = 8, height = 12)

for (state_key in names(dge_type3_results)) {
  state <- gsub("_", " ", state_key)
  
  # Get cells for this state
  state_cells <- Cells(pdos_matched)[pdos_matched$state == state | 
                                      pdos_matched$state == state_key]
  if (length(state_cells) < 100) next
  
  state_obj <- subset(pdos_matched, cells = state_cells)
  
  ht <- make_type3_heatmap(dge_type3_results[[state_key]], state_obj, state)
  if (!is.null(ht)) draw(ht)
}

dev.off()

####################
# TYPE 3: Enrichment analysis (Hallmark + GO)
####################

message("\n=== TYPE 3: Enrichment analysis per state ===")

enrich_type3_up <- list()
enrich_type3_down <- list()

for (state_key in names(dge_type3_up)) {
  message(sprintf("  Enrichment for %s (Up)", state_key))
  enrich_type3_up[[state_key]] <- run_enrichment(dge_type3_up[[state_key]], paste0("Type3_", state_key, "_up"))
  
  message(sprintf("  Enrichment for %s (Down)", state_key))
  enrich_type3_down[[state_key]] <- run_enrichment(dge_type3_down[[state_key]], paste0("Type3_", state_key, "_down"))
}

# Save enrichment results
saveRDS(list(up = enrich_type3_up, down = enrich_type3_down),
        file.path(out_dir, "DGE_Type3_enrichment.rds"))

# Enrichment heatmaps
pdf(file.path(out_dir, "Type3_enrichment_heatmaps.pdf"), width = 12, height = 10)

if (length(enrich_type3_up) > 0) {
  plot_enrich_heatmap(enrich_type3_up, "Hallmark", "Type 3 Upregulated (All patients per state)")
  plot_enrich_heatmap(enrich_type3_up, "GO", "Type 3 Upregulated (All patients per state)")
}

if (length(enrich_type3_down) > 0) {
  plot_enrich_heatmap(enrich_type3_down, "Hallmark", "Type 3 Downregulated (All patients per state)")
  plot_enrich_heatmap(enrich_type3_down, "GO", "Type 3 Downregulated (All patients per state)")
}

dev.off()

####################
# TYPE 3: Summary statistics
####################

message("\n=== TYPE 3: Generating summary ===")

type3_summary <- data.frame(
  State = gsub("_", " ", names(dge_type3_results)),
  Total_DEGs = sapply(dge_type3_results, function(df) {
    sum(df$p_val_adj < padj_threshold & abs(df$avg_log2FC) > logfc_threshold)
  }),
  Upregulated = sapply(dge_type3_up[names(dge_type3_results)], length),
  Downregulated = sapply(dge_type3_down[names(dge_type3_results)], length)
)

write.csv(type3_summary, file.path(out_dir, "Type3_summary.csv"), row.names = FALSE)

# Type 3 summary barplot
pdf(file.path(out_dir, "Type3_summary_barplot.pdf"), width = 10, height = 6)

p3 <- type3_summary %>%
  pivot_longer(cols = c(Upregulated, Downregulated), names_to = "Direction", values_to = "Count") %>%
  ggplot(aes(State, Count, fill = Direction)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c(Upregulated = "firebrick3", Downregulated = "steelblue")) +
  labs(title = "Type 3: DEG counts per state (All patients combined)",
       x = NULL, y = "Number of DEGs") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
  )

print(p3)

dev.off()

message("\n=== DONE ===")
message(sprintf("Output directory: %s", file.path(getwd(), out_dir)))
