library(gficf)
library(ggplot2)

merged_obj <- readRDS("PDOs_merged.rds")

data <- gficf(M = merged_obj@assays$RNA$counts)
data <- runPCA(data = data,dim = 10,use.odgenes = T)
data <- runReduction(data = data,reduction = "umap",nt = 2,verbose = T,n_neighbors=150)

data$embedded$ccl = sapply(
  strsplit(x = rownames(data$embedded),
           split = "_",fixed = T)
  ,function(x) x[1]
)

p = plotCells(data = data,pointShape = 19,colorBy = "ccl") + 
  xlab("UMAP 1") + 
  ylab("UMAP 2")

plot(p)

data <- runScGSEA(
  data = data,
  geneID = "symbol",  # or "symbol" if you're using gene symbols
  species = "human", 
  pathway.list = gene_programs, 
  nmf.k = 100,
  fdr.th = 0.1,
  rescale = "none",
  verbose = TRUE
)


data = clustcellsBYscGSEA(data,
                          method = "fgraph",
                          pca = 10,
                          k = 10,
                          resolution = .05,
                          n.start = 10,
                          n.iter = 50,
                          verbose = T)

p = plotCells(data = data,pointShape = 19,colorBy = "cluster.by.scGSEA") + 
  xlab("UMAP 1") + 
  ylab("UMAP 2")

plot(p)



score_matrix <- t(as.matrix(data$scgsea$x))

sample_ids <- unique(data$embedded$ccl)
pseudobulk_list <- sapply(sample_ids, function(id) {
  cells_in_sample <- rownames(data$embedded[data$embedded$ccl == id, ])
  rowMeans(score_matrix[, cells_in_sample], na.rm = TRUE)
}, simplify = FALSE, USE.NAMES = TRUE) # USE.NAMES keeps the sample IDs as names
score_matrix <- do.call(cbind, pseudobulk_list)




library(Seurat)

genes_to_keep <- all_genes[all_genes %in% rownames(merged_obj)]
subset_obj <- subset(merged_obj, features = genes_to_keep)

subset_obj <- ScaleData(subset_obj, verbose = FALSE)
subset_obj <- RunPCA(subset_obj, features = genes_to_keep, verbose = FALSE)
subset_obj <- RunUMAP(subset_obj, dims = 1:30) # Adjust dims as needed

DimPlot(subset_obj, reduction = "umap")





module_names <- c(
  "Mitosis1", "DNA repair2", "Chromatin Assembly3", "Muscle Hypertrophy4",
  "NA...(SUR735)5", "Epithelial Differentiation6", "Wnt Signaling7", "NA...(SUR759)8",
  "Epithelial Proliferation9", "DNA Regulation10"
)
module_scores <- subset_obj@meta.data[, module_names]
module_assay <- CreateAssayObject(counts = t(module_scores))
module_obj <- CreateSeuratObject(counts = module_assay, meta.data = subset_obj@meta.data)

module_obj <- ScaleData(module_obj)
module_obj <- RunPCA(module_obj, features = rownames(module_obj), approx = FALSE) # Use all 10 features for PCA
module_obj <- RunUMAP(module_obj, dims = 1:10) # Use all 10 PCs

p1 <- DimPlot(module_obj, reduction = "umap", group.by = "orig.ident", label = FALSE)
p2 <- DimPlot(module_obj, reduction = "umap", group.by = "Reponse_mandard", label = FALSE)
p3 <- DimPlot(module_obj, reduction = "umap", group.by = "Tumour_type", label = FALSE)
combined_plot <- p1 + p2 + p3 + plot_layout(ncol = 3)
combined_plot

cell_level_metadata <- annotation_df[match(module_obj$orig.ident, rownames(annotation_df)), ]

module_obj <- AddMetaData(
  object = module_obj,
  metadata = cell_level_metadata
)
colnames(module_obj@meta.data)[colnames(module_obj@meta.data) == "Response based on Mandard"] <- "Reponse_mandard"
colnames(module_obj@meta.data)[colnames(module_obj@meta.data) == "Tumour type: Oesophageal, GOJ type I-III, gastric"] <- "Tumour_type"

