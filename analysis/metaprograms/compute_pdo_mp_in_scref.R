####################
# Score PDO MPs in scRef cells using UCell
#
# Input:
#   PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds
#   scRef_Pipeline/ref_outs/EAC_Ref_epi.rds
#
# Output:
#   scRef_Pipeline/ref_outs/UCell_pdo_in_scref.rds
####################

library(Seurat)
library(UCell)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# Load PDO metaprograms
####################
message("Loading PDO metaprograms...")
MP_pdo <- readRDS("Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds")

# Get MP genes and filter bad MPs (silhouette < 0)
mp.genes <- MP_pdo$metaprograms.genes
bad_mps <- which(MP_pdo$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
coverage_tbl <- MP_pdo$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]

# Filter MPs
mp.genes <- mp.genes[!names(mp.genes) %in% c(bad_mp_names, low_coverage_mps)]
message("Retained PDO MPs: ", paste(names(mp.genes), collapse = ", "))

####################
# Load scRef Seurat object
####################
message("Loading scRef Seurat object...")
scref <- readRDS("/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds")
message("scRef cells: ", ncol(scref))

####################
# Score PDO MPs in scRef cells
####################
message("Computing UCell scores for PDO MPs in scRef cells...")
scref <- AddModuleScore_UCell(scref, features = mp.genes, ncores = 1, name = "")

# Extract MP scores
mp_cols <- grep("^MP", colnames(scref@meta.data), value = TRUE)
ucell_scores <- t(scref@meta.data[, mp_cols, drop = FALSE])

message("UCell scores dimensions: ", dim(ucell_scores))

####################
# Save
####################
saveRDS(ucell_scores, file = "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_pdo_in_scref.rds")

message("=== DONE ===")
message("Saved: /rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/UCell_pdo_in_scref.rds")
