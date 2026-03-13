####################
# Auto_update_optimal_mp.R
# Re-run post-NMF analysis with the optimal nMP determined by inflection analysis
# Updates: MP_outs_default.rds, GO_outs.rds, PDOs_final.rds, VlnPlots
# Run with gnmf env (needs GeneNMF, UCell, Seurat)
####################

library(GeneNMF)
library(RColorBrewer)
library(UCell)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# Set optimal nMP (from Auto_find_optimal_nmf.R inflection analysis)
####################
optimal_nMP <- 13

message(paste0("Using optimal nMP = ", optimal_nMP))

# Load optimal metaprogram result
geneNMF.metaprograms <- readRDS(file.path("Metaprogrammes_Results",
                                          paste0("geneNMF_metaprograms_nMP_", optimal_nMP, ".rds")))
saveRDS(geneNMF.metaprograms, file = "MP_outs_default.rds")
message("Saved MP_outs_default.rds")

# Heatmap for optimal nMP
n_colors <- min(optimal_nMP, 12)
anno_colors <- brewer.pal(n = max(n_colors, 3), name = "Paired")
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
message("Saved metaprograms_heatmap.png")

####################
# GSEA per metaprogram (GO:BP)
####################
pdos <- readRDS("PDOs_merged.rds")
pdos <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")

mp.genes <- geneNMF.metaprograms$metaprograms.genes

top_p <- lapply(mp.genes, function(program) {
  runGSEA(program, universe = rownames(pdos), category = "C5", subcategory = "GO:BP")
})
saveRDS(top_p, "GO_outs.rds")
message("Saved GO_outs.rds")

####################
# UCell scoring
####################
pdos <- AddModuleScore_UCell(pdos, features = mp.genes, ncores = 4, name = "")
saveRDS(pdos, "PDOs_final.rds")
message("Saved PDOs_final.rds")

####################
# Violin plots
####################
png("vln_origident.png",
    width = 5000, height = 2500, res = 300)
VlnPlot(
  pdos,
  features = names(mp.genes),
  group.by = "orig.ident",
  pt.size = 0,
  ncol = 5
)
dev.off()

png("vln_clinical_response.png",
    width = 3500, height = 2500, res = 300)
VlnPlot(
  pdos,
  features = names(mp.genes),
  group.by = "Clinical.response.at.OG.MDT..responder.non.responder",
  pt.size = 0,
  ncol = 5
)
dev.off()

png("vln_batch.png",
    width = 3500, height = 2500, res = 300)
VlnPlot(
  pdos,
  features = names(mp.genes),
  group.by = "Batch",
  pt.size = 0,
  ncol = 5
)
dev.off()

message("All outputs updated with optimal nMP = ", optimal_nMP)
