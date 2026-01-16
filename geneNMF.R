library(GeneNMF)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(UCell)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")
pdos.list <- readRDS("PDOs_list_PDOs.rds")

geneNMF.programs <- multiNMF(pdos.list, assay="RNA", k=4:9, min.exp = 0.05)
saveRDS(geneNMF.programs, file="geneNMF_outs.rds")

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        specificity.weight = 5,
                                        weight.explained = 0.5,
                                        nMP=10, 
                                        min.confidence = 0.5)

anno_colors <- brewer.pal(n=10, name="Paired")
names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

png("metaprograms_heatmap.png",
    width = 3000,       # pixels
    height = 2500,      # pixels
    res = 300)          # DPI (good for heatmaps)
plotMetaPrograms(
  geneNMF.metaprograms,
  annotation_colors = anno_colors,
  similarity.cutoff = c(0,1)
)
dev.off()

pdos <- readRDS("PDOs_merged.rds")
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(pdos), category = "C5", subcategory = "GO:BP")
})
saveRDS(top_p, "GO_outs.rds")

mp.genes <- geneNMF.metaprograms$metaprograms.genes
pdos <- AddModuleScore_UCell(pdos, features = mp.genes, ncores=4, name = "")

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