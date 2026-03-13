library(GeneNMF)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(UCell)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
pdos.list <- readRDS("PDOs_list_PDOs.rds")
pdos.list$SUR843T3_PDO <- NULL

geneNMF.programs <- multiNMF(pdos.list, assay="RNA", k=4:9, min.exp = 0.05)
saveRDS(geneNMF.programs, file="geneNMF_outs.rds")

####################
# Loop over nMP values (4:20) to find optimal number of metaprograms
# Each result saved to Metaprogrammes_Results/ for downstream silhouette/WSS analysis
####################
dir.create("Metaprogrammes_Results", showWarnings = FALSE)
k_vals <- 4:20
for (k in k_vals) {
  message(paste0("Running getMetaPrograms with nMP = ", k))
  geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                          metric = "cosine",
                                          specificity.weight = 5,
                                          weight.explained = 0.5,
                                          nMP = k,
                                          min.confidence = 0.5)
  saveRDS(geneNMF.metaprograms, file = file.path("Metaprogrammes_Results",
                                                  paste0("geneNMF_metaprograms_nMP_", k, ".rds")))

  # Save heatmap for each nMP
  n_colors <- min(k, 12)
  anno_colors <- brewer.pal(n = max(n_colors, 3), name = "Paired")
  anno_colors <- anno_colors[seq_len(length(geneNMF.metaprograms$metaprograms.genes))]
  names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

  png(file.path("Metaprogrammes_Results", paste0("metaprograms_heatmap_nMP_", k, ".png")),
      width = 3000, height = 2500, res = 300)
  plotMetaPrograms(geneNMF.metaprograms,
                   annotation_colors = anno_colors,
                   similarity.cutoff = c(0, 1))
  dev.off()
  message(paste0("Saved nMP = ", k))
}

####################
# Optimal nMP=13 determined by inflection-point analysis (Auto_find_optimal_nmf.R)
# Silhouette inflection: 13, WSS elbow: 16
####################
default_nMP <- 13
geneNMF.metaprograms <- readRDS(file.path("Metaprogrammes_Results",
                                          paste0("geneNMF_metaprograms_nMP_", default_nMP, ".rds")))
saveRDS(geneNMF.metaprograms, file = "MP_outs_default.rds")

anno_colors <- brewer.pal(n = min(default_nMP, 12), name = "Paired")
anno_colors <- anno_colors[seq_len(length(geneNMF.metaprograms$metaprograms.genes))]
names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

png("metaprograms_heatmap.png",
    width = 3000, height = 2500, res = 300)
plotMetaPrograms(
  geneNMF.metaprograms,
  annotation_colors = anno_colors,
  similarity.cutoff = c(0, 1)
)
dev.off()
####################

pdos <- readRDS("PDOs_merged.rds")
pdos <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(pdos), category = "C5", subcategory = "GO:BP")
})
saveRDS(top_p, "GO_outs.rds")

mp.genes <- geneNMF.metaprograms$metaprograms.genes
pdos <- AddModuleScore_UCell(pdos, features = mp.genes, ncores=4, name = "")
saveRDS(pdos, "PDOs_final.rds")

png("vln_origident.png",
    width = 5000,   # wide
    height = 2500,
    res = 300)
VlnPlot(
  pdos,
  features = names(mp.genes),
  group.by = "orig.ident",
  pt.size = 0,
  ncol = 5
)
dev.off()

png("vln_clinical_response.png",
    width = 3500,
    height = 2500,
    res = 300)
VlnPlot(
  pdos,
  features = names(mp.genes),
  group.by = "Clinical.response.at.OG.MDT..responder.non.responder",
  pt.size = 0,
  ncol = 5
)
dev.off()

png("vln_batch.png",
    width = 3500,
    height = 2500,
    res = 300)
VlnPlot(
  pdos,
  features = names(mp.genes),
  group.by = "Batch",
  pt.size = 0,
  ncol = 5
)
dev.off()
