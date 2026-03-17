####################
# Auto_ucell_vlnplot.R
# UCell scoring and violin plots for optimal nMP
# Separate from Auto_update_optimal_mp.R because UCell is memory-intensive
# Run with gnmf env via PBS
####################

library(UCell)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# Load optimal MP result and merged Seurat
####################
geneNMF.metaprograms <- readRDS("MP_outs_default.rds")
mp.genes <- geneNMF.metaprograms$metaprograms.genes

message(paste0("Loaded ", length(mp.genes), " metaprograms"))
message(paste0("MP names: ", paste(names(mp.genes), collapse = ", ")))

pdos <- readRDS("PDOs_merged.rds")
pdos <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")
message(paste0("Loaded ", ncol(pdos), " cells"))

####################
# UCell scoring
####################
message("Starting UCell scoring...")
pdos <- AddModuleScore_UCell(pdos, features = mp.genes, ncores = 4, name = "")
message("UCell scoring complete")

# Save UCell scores separately for lightweight downstream use
ucell_cols <- grep("^MP", colnames(pdos@meta.data), value = TRUE)
ucell_scores <- pdos@meta.data[, ucell_cols, drop = FALSE]
saveRDS(ucell_scores, "UCell_scores_filtered.rds")
message("Saved UCell_scores_filtered.rds")

saveRDS(pdos, "PDOs_final.rds")
message("Saved PDOs_final.rds")

####################
# Violin plots
####################
message("Generating violin plots...")
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
message("Saved vln_origident.png")

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
message("Saved vln_clinical_response.png")

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
message("Saved vln_batch.png")

message("All UCell + VlnPlot outputs complete")
