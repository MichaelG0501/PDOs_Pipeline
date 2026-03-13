library(pheatmap)
library(viridis)

sim_matrix <- geneNMF.metaprograms$programs.similarity
nmf_tree <- geneNMF.metaprograms$programs.tree
mp_clusters <- geneNMF.metaprograms$programs.clusters
sample_names <- sub("\\..*$", "", rownames(sim_matrix))
annotation_df <- data.frame(
  Sample = sample_names,
  row.names = rownames(sim_matrix)
)
annotation_df$Batch <- ifelse(
  grepl("_PDO", annotation_df$Sample) & 
    !grepl("Untreated_PDO|Treated_PDO", annotation_df$Sample),
  "batch_cynthia",
  ifelse(grepl("Untreated_PDO|Treated_PDO", annotation_df$Sample), "batch2", NA)
)
annotation_df$Metaprogram <- as.factor(mp_clusters)
sample_levels <- unique(annotation_df$Sample)
sample_palette <- setNames(viridis(length(sample_levels), option = "D"), sample_levels)
mp_levels <- levels(annotation_df$Metaprogram)
mp_palette <- setNames(
  RColorBrewer::brewer.pal(n = max(3, length(mp_levels)), name = "Paired")[1:length(mp_levels)],
  mp_levels
)

original_rnames <- rownames(annotation_df)
PDOs_all_meta <- readRDS("PDOs_all_meta.rds")
resp_lookup <- PDOs_all_meta %>%
  dplyr::select(orig.ident, Response.based.on.Mandard) %>%
  distinct()
annotation_df <- annotation_df %>%
  left_join(resp_lookup, by = c("Sample" = "orig.ident"))
rownames(annotation_df) <- original_rnames
annotation_df$Response.based.on.Mandard <- dplyr::recode(
  annotation_df$Response.based.on.Mandard,
  "R" = "Responder",
  "NR" = "Non-Responder"
)
annotation_df$Response.based.on.Mandard[is.na(annotation_df$Response.based.on.Mandard)] <- "NA"
response_palette <- c("Responder" = "darkorange", 
                      "Non-Responder" = "navy", 
                      "NA" = "grey90")
anno_colors <- list(
  Batch = c(batch_cynthia = "brown", batch2 = "darkgreen"),
  Sample = sample_palette,
  Metaprogram = mp_palette, 
  Response.based.on.Mandard = response_palette
)

annotation_df <- annotation_df[, c("Metaprogram", "Batch", "Response.based.on.Mandard"), drop = FALSE]
pheatmap(
  sim_matrix,
  cluster_rows = nmf_tree,
  cluster_cols = nmf_tree,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  annotation_colors = anno_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Metaprogram Similarity (GeneNMF Tree Order)",
  color = viridis(100, option = "A", direction = -1), 
  fontsize = 8,           # Smaller global font size
  fontsize_number = 7,
)
#########################
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(viridis)
library(RColorBrewer)

# ---- 1. Data Preparation & Cleaning ----
sim_matrix  <- geneNMF.metaprograms$programs.similarity
nmf_tree    <- geneNMF.metaprograms$programs.tree
mp_clusters <- geneNMF.metaprograms$programs.clusters

# Prepare Annotation Data
sample_names <- sub("\\..*$", "", rownames(sim_matrix))
annotation_df <- data.frame(
  Sample = sample_names, 
  # Add "MP" prefix here
  Metaprogram = paste0("MP", mp_clusters),
  row.names = rownames(sim_matrix),
  stringsAsFactors = FALSE
)

annotation_df$Study <- sapply(strsplit(annotation_df$Sample, "_"), function(x) paste(x[1], x[2], sep="_"))
rownames(annotation_df) <- rownames(sim_matrix)

# 2. FILTER: Remove NA, "MPNA", or empty
# Note: Since we pasted "MP", NA becomes "MPNA"
keep_names <- rownames(annotation_df)[!is.na(mp_clusters) & mp_clusters != ""]

# 3. Apply Tree Order
# Get names in tree order, then keep only those in our 'keep' list
all_ordered_names <- nmf_tree$labels[nmf_tree$order]
final_ordered_names <- all_ordered_names[all_ordered_names %in% keep_names]

# Subset and Reorder
sim_matrix_ord <- sim_matrix[final_ordered_names, final_ordered_names]
annotation_ord <- annotation_df[final_ordered_names, , drop = FALSE]

# Lock factor levels to the tree sequence to prevent slice shuffling
annotation_ord$Metaprogram <- factor(annotation_ord$Metaprogram, levels = unique(annotation_ord$Metaprogram))

# ---- 4. Colors & Aesthetics ----
mp_palette <- setNames(
  colorRampPalette(brewer.pal(8, "Paired"))(length(levels(annotation_ord$Metaprogram))),
  levels(annotation_ord$Metaprogram)
)

study_levels <- unique(annotation_ord$study)
study_palette <- setNames(viridis::viridis(length(study_levels), option = "turbo"), study_levels)

anno_colors <- list(Metaprogram = mp_palette, study = study_palette)
col_fun <- colorRamp2(
  c(0.00, 0.12, 0.22, 0.70, 1.00),
  c("#FFFFFF", "#F6E8A6", "#E76F51", "#5E2A84", "#000000")
)
# ---- 5. Construct Heatmap ----
top_ha <- HeatmapAnnotation(
  df = annotation_ord %>% dplyr::select(Metaprogram, study),
  col = anno_colors,
  show_annotation_name = FALSE, # No label for the annotation tracks
  show_legend = TRUE,           # Keeps the legend on the side
  simple_anno_size = unit(3, "mm")
)

left_ha <- rowAnnotation(
  df = annotation_ord %>% dplyr::select(Metaprogram, study),
  col = anno_colors,
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- Heatmap(
  sim_matrix_ord,
  name = "Similarity",
  col = col_fun,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # Split by Metaprogram
  row_split = annotation_ord$Metaprogram,
  column_split = annotation_ord$Metaprogram,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  
  # Visuals: No block borders, small gaps
  rect_gp = gpar(col = NA), 
  border = FALSE,
  row_gap = unit(1, "mm"),
  column_gap = unit(1, "mm"),
  
  # Remove all labels/titles
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_title = NULL,
  column_title = NULL,
  
  top_annotation = top_ha,
  left_annotation = left_ha, 
  use_raster = FALSE
)

# ---- 6. Draw ----
draw(ht, merge_legends = TRUE)

################

module_scores <- readRDS("UCell_default.rds")
pdos@meta.data[, grepl("^MP", colnames(pdos@meta.data))] <- NULL
pdos@meta.data <- cbind(pdos@meta.data, module_scores)
mp.genes <- geneNMF.metaprograms$metaprograms.genes
# update <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")
# matrix <- update@meta.data[,names(mp.genes)]
# dimred <- as.matrix(matrix)
# colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
# update@reductions[["MPsignatures"]] <- new("DimReduc",
#                                             cell.embeddings = dimred,
#                                             assay.used = "RNA",
#                                             key = "MP_",
#                                             global = FALSE)
# update <- RunUMAP(update, reduction="MPsignatures", dims=1:length(update@reductions[["MPsignatures"]]),
#                    metric = "euclidean", reduction.name = "umap_MP")
# update$Batch <- ifelse(
#   grepl("_PDO", update$orig.ident) & 
#     !grepl("Untreated_PDO|Treated_PDO", update$orig.ident),
#   "batch_cynthia",
#   ifelse(grepl("Untreated_PDO|Treated_PDO", update$orig.ident), "batch2", NA)
# )

cells_to_keep <- rownames(pdos@meta.data)[!is.na(pdos@meta.data$Response.based.on.Mandard)]
subset_pdos <- subset(pdos, cells = cells_to_keep)
subset_pdos$Response.based.on.Mandard <- dplyr::recode(
  subset_pdos$Response.based.on.Mandard,
  "R" = "Responder",
  "NR" = "Non-Responder"
)
subset_pdos$Response.based.on.Mandard <- factor(
  subset_pdos$Response.based.on.Mandard,
  levels = c("Responder", "Non-Responder")
)
VlnPlot(subset_pdos, features=names(mp.genes)[c(3,6,7)], group.by = "Response.based.on.Mandard",
        pt.size = 0, ncol=3)
VlnPlot(subset_pdos, features=names(mp.genes)[c(3,6,7)], group.by = "Mandard.tumour.regression.score",
        pt.size = 0, ncol=3)

subset_pdos$Batch <- ifelse(
  grepl("_PDO", subset_pdos$orig.ident) & 
    !grepl("Untreated_PDO|Treated_PDO", subset_pdos$orig.ident),
  "batch_cynthia",
  ifelse(grepl("Untreated_PDO|Treated_PDO", subset_pdos$orig.ident), "batch2", NA)
)
batch_B <- subset(subset_pdos, subset = Batch == "batch2")
batch_B$Treatment <- ifelse(
  grepl("Treated_PDO", batch_B$orig.ident), "Treated",
  ifelse(grepl("Untreated_PDO", batch_B$orig.ident), "Untreated", NA)
)
VlnPlot(batch_B, features=names(mp.genes)[c(3,6,7)], group.by = "Response.based.on.Mandard",
        pt.size = 0, ncol=3)
VlnPlot(batch_B, features=names(mp.genes)[c(3,6,7)], group.by = "Mandard.tumour.regression.score",
        pt.size = 0, ncol=3)

matrix <- batch_B@meta.data[,names(mp.genes)]
dimred <- as.matrix(matrix)
colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
batch_B@reductions[["MPsignatures"]] <- new("DimReduc",
                                            cell.embeddings = dimred,
                                            assay.used = "RNA",
                                            key = "MP_",
                                            global = FALSE)
batch_B <- RunUMAP(batch_B, reduction="MPsignatures", dims=1:length(batch_B@reductions[["MPsignatures"]]),
                   metric = "euclidean", reduction.name = "umap_MP")
DimPlot(batch_B, reduction = "umap_MP", group.by = "Response.based.on.Mandard") + theme(aspect.ratio = 1)
DimPlot(batch_B, reduction = "umap_MP", group.by = "orig.ident") + theme(aspect.ratio = 1)


batch_A <- subset(subset_pdos, subset = Batch == "batch_cynthia")
batch_A <- subset(batch_A, subset = orig.ident != "SUR843T3_PDO")
VlnPlot(batch_A, features=names(mp.genes)[c(3,6,7)], group.by = "Response.based.on.Mandard",
        pt.size = 0, ncol=3)
VlnPlot(batch_A, features=names(mp.genes)[c(3,6,7)], group.by = "Mandard.tumour.regression.score",
        pt.size = 0, ncol=3)

matrix <- batch_A@meta.data[,names(mp.genes)]
dimred <- as.matrix(matrix)
colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
batch_A@reductions[["MPsignatures"]] <- new("DimReduc",
                                        cell.embeddings = dimred,
                                        assay.used = "RNA",
                                        key = "MP_",
                                        global = FALSE)
batch_A <- RunUMAP(batch_A, reduction="MPsignatures", dims=1:length(batch_A@reductions[["MPsignatures"]]),
               metric = "euclidean", reduction.name = "umap_MP")
DimPlot(batch_A, reduction = "umap_MP", group.by = "Response.based.on.Mandard") + theme(aspect.ratio = 1)
DimPlot(batch_A, reduction = "umap_MP", group.by = "orig.ident") + theme(aspect.ratio = 1)


#########################################

pdos <- subset(pdos, subset = Batch != "PDO")
pdos <- subset(pdos, subset = SUR != "SUR1121" & SUR != "SUR1141")


pdos <- FindVariableFeatures(pdos)
pdos <- ScaleData(pdos)
pdos <- RunPCA(pdos)
pdos <- FindNeighbors(pdos, dims = 1:50)

pdos <- FindClusters(pdos, resolution = 0.8, algorithm = 1)
pdos$leiden_clusters <- Idents(pdos)

pdos <- RunUMAP(pdos, dims = 1:50)

library(readxl)
library(dplyr)
data <- read_excel("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx", sheet = 1)
data <- as.matrix(data)
rownames(data) <- data[, 1]
data <- t(data[ , -1])
rownames(data) <- paste0("SUR", rownames(data))
data <- data.frame(SUR = rownames(data), data, row.names = NULL)
cell_names <- colnames(pdos)
new_cols <- setdiff(colnames(data), "SUR")
pdos@meta.data <- pdos@meta.data %>%
  select(-any_of(new_cols)) %>%   # remove old versions
  left_join(data, by = "SUR")
rownames(pdos@meta.data) <- cell_names

p1 <- DimPlot(pdos, group.by = "orig.ident", label = TRUE) + ggtitle("Louvain Clustering")
p2 <- DimPlot(pdos, group.by = "Collection timepoint", label = FALSE) + ggtitle("orig.ident")
p3 <- DimPlot(pdos, group.by = "Batch", label = FALSE) + ggtitle("Batch")
p4 <- DimPlot(pdos, group.by = "Clinical.response.at.OG.MDT..responder.non.responder", label = FALSE) + ggtitle("Response")

combined_plot <- (p1 | p2) / (p3 | p4)

######################################

mp.genes <- geneNMF.metaprograms$metaprograms.genes
pdos$`Clinical.response.at.OG.MDT..responder.non.responder` <- factor(
  pdos$`Clinical.response.at.OG.MDT..responder.non.responder`,
  levels = c("Responder", "Non-responder")
)
pdos$Batch <- factor(
  pdos$Batch,
  levels = c("Untreated_PDO", "Treated_PDO")
)
VlnPlot(pdos, features=names(mp.genes)[c(3,6,7)], group.by = "Clinical.response.at.OG.MDT..responder.non.responder",
        pt.size = 0, ncol=3)
VlnPlot(pdos, features=names(mp.genes)[c(3,6,7)], group.by = "Collection timepoint",
        pt.size = 0, ncol=3)
VlnPlot(pdos, features=names(mp.genes)[c(3,6,7)], group.by = "Batch",
        pt.size = 0, ncol=3)

matrix <- pdos@meta.data[,names(mp.genes)]
dimred <- as.matrix(matrix)
colnames(dimred) <- paste0("MP_",seq(1, ncol(dimred)))
pdos@reductions[["MPsignatures"]] <- new("DimReduc",
                                         cell.embeddings = dimred,
                                         assay.used = "RNA",
                                         key = "MP_",
                                         global = FALSE)
pdos <- RunUMAP(pdos, reduction="MPsignatures", dims=1:length(pdos@reductions[["MPsignatures"]]),
                metric = "euclidean", reduction.name = "umap_MP")
DimPlot(pdos, reduction = "umap_MP", group.by = "Clinical.response.at.OG.MDT..responder.non.responder") + theme(aspect.ratio = 1)
DimPlot(pdos, reduction = "umap_MP", group.by = "orig.ident") + theme(aspect.ratio = 1)
DimPlot(pdos, reduction = "umap_MP", group.by = "Batch") + theme(aspect.ratio = 1)

FeaturePlot(pdos, features = names(mp.genes)[3])

library(viridis)
FeaturePlot(pdos, features = names(mp.genes), reduction = "umap_MP", ncol=4) &
  scale_color_viridis(option="B") &
  theme(aspect.ratio = 1, axis.text=element_blank(), axis.ticks=element_blank())



library(viridis)
library(patchwork)

plot_list <- FeaturePlot(
  pdos,
  features = names(mp.genes),
  combine = FALSE 
) 

plot_list <- lapply(plot_list, function(x) {
  x + scale_color_viridis(option = "B") +
    theme(aspect.ratio = 1, 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title = element_blank(), # Remove axis titles to save space
          plot.title = element_text(size = 7, margin = margin(b = 1)), # Much smaller title
          plot.margin = margin(1, 1, 1, 1), # Minimal padding
          legend.position = "none")
})

# 2. DimPlots: Match the styling
p_orig <- DimPlot(pdos, group.by = "Clinical.response.at.OG.MDT..responder.non.responder") + 
  ggtitle("Response") + 
  theme(aspect.ratio = 1, 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 7, margin = margin(b = 1)),
        plot.margin = margin(1, 1, 1, 1),
        legend.position = "none")

p_batch <- DimPlot(pdos, group.by = "Batch") + 
  ggtitle("Treatment status") + 
  theme(aspect.ratio = 1, 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 7, margin = margin(b = 1)),
        plot.margin = margin(1, 1, 1, 1),
        legend.position = "none")

# 3. Combine into tight 3x4 grid
wrap_plots(c(plot_list, list(p_orig, p_batch)), ncol = 4) & 
  theme(plot.margin = margin(0, 0, 0, 0))