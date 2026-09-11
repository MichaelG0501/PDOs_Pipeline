####################
# Analysis registry:
#   Status: active upstream score-regeneration utility
#   Script: analysis/metaprograms/centred/Auto_07a_ucell_scoring.R
#   Methodology: analysis/methodology/metaprograms/centred/Auto_centred_ordering_and_state_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Recomputes rank-based UCell scores for the final centred-refined PDO MP
#     gene lists and the external 3CA MP signatures on all retained PDO cells.
#
#   Inputs:
#     - live: /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv
#     - live: centred_mp_refinement/merged_refined_mp_genes.rds
#     - live: PDOs_outs/PDOs_merged.rds
#
#   Outputs (live: PDOs_outs/):
#     - centred_mp_refinement/merged_refined_ucell_scores.rds
#     - UCell_3CA_MPs.rds
#   Outputs (ephemeral cache):
#     - centred_mp_refinement/intermediate/merged_refined_ucell_scores.rds
#   Downstream use:
#     Both score matrices are persistent inputs to state definition and current
#     centred correlation/comparison workflows.
#   Cache/replot behavior: always recomputes scores; no cache-read branch.
#   Run command: qsub Auto_run_centred_07a_ucell.sh
#
#   Conda env: gnmf
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
})

live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
ephemeral_base <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"

mp_genes_file <- file.path(live_base, "centred_mp_refinement", "merged_refined_mp_genes.rds")
if (!file.exists(mp_genes_file)) {
  mp_genes_file <- file.path(ephemeral_base, "centred_mp_refinement", "intermediate", "merged_refined_mp_genes.rds")
}
nmfs_csv_file <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
seurat_file <- file.path(live_base, "PDOs_merged.rds")

cat("Loading merged refined MP genes...\n")
refined_genes <- readRDS(mp_genes_file)

cat("Loading 3CA MP genes...\n")
nmfs_csv <- read.csv(nmfs_csv_file, stringsAsFactors = FALSE)
nmfs_list <- lapply(as.list(nmfs_csv), function(x) unique(x[x != "" & !is.na(x)]))

clean_3ca_name <- function(x) {
  x <- gsub("^X?MP", "3CA_MP", x)
  x <- gsub("\\.", " ", x)
  x
}
names(nmfs_list) <- clean_3ca_name(names(nmfs_list))

cat("Loading PDOs merged Seurat object...\n")
pdos_merged <- readRDS(seurat_file)

all_features <- c(refined_genes, nmfs_list)

cat("Running ScoreSignatures_UCell...\n")
counts_matrix <- GetAssayData(pdos_merged, layer = "data")
ucell_matrix <- ScoreSignatures_UCell(counts_matrix, features = all_features, ncores = 4)

# ScoreSignatures_UCell appends "_UCell" to colnames, let's clean column names
colnames(ucell_matrix) <- sub("_UCell$", "", colnames(ucell_matrix))

ucell_df <- as.data.frame(ucell_matrix)

ucell_refined <- ucell_df[, names(refined_genes), drop = FALSE]
ucell_3ca <- ucell_df[, names(nmfs_list), drop = FALSE]

outdir_refined <- file.path(live_base, "centred_mp_refinement")
dir.create(outdir_refined, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ephemeral_base, "centred_mp_refinement", "intermediate"), recursive = TRUE, showWarnings = FALSE)

saveRDS(ucell_refined, file.path(outdir_refined, "merged_refined_ucell_scores.rds"))
saveRDS(ucell_refined, file.path(ephemeral_base, "centred_mp_refinement", "intermediate", "merged_refined_ucell_scores.rds"))

saveRDS(ucell_3ca, file.path(live_base, "UCell_3CA_MPs.rds"))

cat("Saved UCell scores cleanly to live and ephemeral.\n")
