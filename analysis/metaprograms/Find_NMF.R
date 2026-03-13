library(Seurat)
library(NMF)
library(reshape2)
library(ggplot2)
library(scales)
library(doParallel)
library(foreach)


n_cores <- 11
cl <- makeCluster(n_cores)
registerDoParallel(cl)

nmf.options(parallel = 5)

tmdata_list <- readRDS("PDOs_all.rds")

clusterExport(cl, c("tmdata_list"))

Genes_nmf_w_basis_list <- foreach(name = names(tmdata_list),
                                  .packages = c("Seurat", "NMF", "ggplot2", "scales"),
                                  .errorhandling = 'pass') %dopar% {
                                    tryCatch({
                                      obj <- tmdata_list[[name]]
                                      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
                                      var_genes <- VariableFeatures(obj)
                                      
                                      expr_matrix <- GetAssayData(obj, layer = "CPM")
                                      expr_var <- expr_matrix[var_genes, ] / 10
                                      expr_norm <- t(scale(t(expr_var), center = TRUE, scale = FALSE))
                                      expr_norm <- pmax(expr_norm, 0)
                                      
                                      rank_all <- list()
                                      for (rank in 4:9) {
                                        message(paste("Running NMF for", name, "rank", rank))
                                        nmf_result <- nmf(expr_norm, rank = rank, nrun = 10, method = "brunet")
                                        basis_matrix <- basis(nmf_result)
                                        colnames(basis_matrix) <- paste0(name, "_rank4_9_nrun10.RDS.", rank, ".", seq_len(rank))
                                        rank_all[[rank]] <- basis_matrix
                                        gc()
                                      }
                                      
                                      sample_name <- paste0(name, "_rank4_9_nrun10.RDS")
                                      result <- list()
                                      result[[sample_name]] <- do.call(cbind, rank_all)
                                      result
                                      
                                    }, error = function(e) {
                                      warning(paste("Error in sample:", name, "->", e$message))
                                      NULL
                                    })
                                  }

Genes_nmf_w_basis <- do.call(c, Genes_nmf_w_basis_list)

saveRDS(Genes_nmf_w_basis, "Genes_nmf_w_basis.rds")
stopCluster(cl)

Genes_nmf_w_basis <- readRDS("Genes_nmf_w_basis.rds")




robust_nmf_programs <- function(nmf_programs, intra_min = 35, intra_max = 10) {
  
  # Select NMF programs based on the minimum overlap with other NMF programs from the same cell line
  intra_intersect <- lapply(nmf_programs, function(z) apply(z, 2, function(x) apply(z, 2, function(y) length(intersect(x,y))))) 
  intra_intersect_max <- lapply(intra_intersect, function(x) apply(x, 2, function(y) sort(y, decreasing = T)[2]))             
  nmf_sel <- lapply(names(nmf_programs), function(x) nmf_programs[[x]][,intra_intersect_max[[x]]>=intra_min]) 
  names(nmf_sel) <- names(nmf_programs)
  
  # Select NMF programs based on i) the maximum overlap with other NMF programs from the same cell line and
  # ii) the minimum overlap with programs from another cell line
  nmf_sel_unlist <- do.call(cbind, nmf_sel)
  inter_intersect <- apply(nmf_sel_unlist , 2, function(x) apply(nmf_sel_unlist , 2, function(y) length(intersect(x,y)))) ## calculating intersection between all programs
  
  final_filter <- NULL 
  a <- inter_intersect
  b <- sort(apply(a, 2, function(x) sort(x, decreasing = TRUE)[2]), decreasing = T) # for each cell line, ranks programs based on their maximum overlap with programs of other cell lines
  if(length(b) > 1) {
    c <- names(b[1]) 
    for(y in 2:length(b)) {
      if(max(inter_intersect[c,names(b[y])]) <= intra_max) c <- c(c,names(b[y])) # selects programs iteratively from top-down. Only selects programs that have a intersection smaller than 10 with a previously selected programs
    }
    final_filter <- c(final_filter, c)
  } else {
    final_filter <- c(final_filter, names(b))
  }
  return(final_filter)                                                      
}

intra_min_parameter <- 30
intra_max_parameter <- 20

names(Genes_nmf_w_basis) <- gsub("nruns", "nrun", names(Genes_nmf_w_basis))

nmf_programs          <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs          <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 
inter_intersect         <- apply(do.call(cbind, nmf_programs) , 2, function(x) apply(do.call(cbind, nmf_programs) , 2, function(y) length(intersect(x,y)))) 

# not filtering inter sample
nmf_filter_ccle       <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter)  
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)

nmf_intersect         <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 

nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

library(pheatmap)

highlight <- colnames(inter_intersect) %in% nmf_filter_ccle
annotation_col <- data.frame(
  Highlight = ifelse(highlight, "Yes", "No"),
  Sample = sub("_rank.*", "", colnames(inter_intersect))
)

annotation_row <- data.frame(
  Highlight = ifelse(rownames(inter_intersect) %in% nmf_filter_ccle, "Yes", "No"),
  Sample = sub("_rank.*", "", rownames(inter_intersect))
)

# Set row names for annotations
rownames(annotation_col) <- colnames(inter_intersect)
rownames(annotation_row) <- rownames(inter_intersect)

# Plot with annotations
pheatmap(
  inter_intersect,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(Highlight = c(Yes = "red", No = "white")),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Highlighted NMF Programs"
)

heatmap(
  nmf_intersect,
  Rowv = nmf_intersect_hc,
  Colv = nmf_intersect_hc,
  scale = "none",
  col = colorRampPalette(c("white", "blue"))(100),
  margins = c(10, 10),
  main = "NMF Program Similarity"
)




overlap_list <- lapply(nmf_filter_ccle, function(prog) {
  overlaps <- inter_intersect[prog, ]
  others <- names(overlaps[overlaps > 20])
  return(others)
})

names(overlap_list) <- nmf_filter_ccle


get_prefix <- function(x) sub("_rank.*", "", x)

unique_sample_counts <- sapply(overlap_list, function(progs) {
  length(unique(get_prefix(progs)))
})
final_nmf <- unique_sample_counts[unique_sample_counts >= 2]

library(clusterProfiler)
library(org.Hs.eg.db)

enrich_res <- lapply(as.data.frame(nmf_programs), function(genes) {
  enrichGO(gene          = genes,
           OrgDb         = org.Hs.eg.db,
           keyType       = "SYMBOL",
           ont           = "BP",
           pAdjustMethod = "BH",
           qvalueCutoff  = 0.05)
})

saveRDS(enrich_res, "enrichGO.rds")

enrich_res_sub <- enrich_res[names(final_nmf)]

top_terms <- lapply(enrich_res_sub, function(df) {
  df <- df[order(df$p.adjust), ]
  head(df$Description, 5)
})

all_terms <- unique(unlist(top_terms))
mat <- sapply(top_terms, function(terms) all_terms %in% terms)
mat <- t(mat) * 1  # binary matrix
rownames(mat) <- names(top_terms)
colnames(mat) <- all_terms

dist_mat <- dist(mat, method = "binary")  # Jaccard-like
hc <- hclust(dist_mat, method = "average")


#################

selected <- unique(unlist(overlap_list[names(final_nmf)], use.names = FALSE))
inter_intersect <- inter_intersect[selected, selected]

highlight <- colnames(inter_intersect) %in% names(final_nmf)
annotation_col <- data.frame(
  Highlight = ifelse(highlight, "Yes", "No"),
  Sample = sub("_rank.*", "", colnames(inter_intersect))
)

annotation_row <- data.frame(
  Highlight = ifelse(rownames(inter_intersect) %in% nmf_filter_ccle, "Yes", "No"),
  Sample = sub("_rank.*", "", rownames(inter_intersect))
)

# Set row names for annotations
rownames(annotation_col) <- colnames(inter_intersect)
rownames(annotation_row) <- rownames(inter_intersect)

# Plot with annotations
pheatmap(
  inter_intersect,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = list(Highlight = c(Yes = "red", No = "white")),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Highlighted NMF Programs"
)


nmf_intersect <- nmf_intersect[names(final_nmf), names(final_nmf)]
nmf_intersect_hc <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

# Compute distance and clustering separately for rows and columns
row_hc <- hclust(as.dist(50 - nmf_intersect), method = "average")
col_hc <- hclust(as.dist(t(50 - nmf_intersect)), method = "average")
row_dend <- reorder(as.dendrogram(row_hc), rowMeans(nmf_intersect))
col_dend <- reorder(as.dendrogram(col_hc), colMeans(nmf_intersect))
nmf_intersect <- nmf_intersect[order.dendrogram(row_dend), order.dendrogram(col_dend)]
heatmap(
  nmf_intersect,
  Rowv = row_dend,
  Colv = col_dend,
  scale = "none",
  col = colorRampPalette(c("white", "blue"))(100),
  margins = c(10, 10),
  main = "NMF Program Similarity"
)



nmf_programs <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
#nmf_programs <- lapply(nmf_programs,toupper) ## convert all genes to uppercase 
nmf_programs <- do.call(cbind, nmf_programs)
all_genes <- nmf_programs[, selected]
all_genes <- as.vector(all_genes)
all_genes <- unique(all_genes)
