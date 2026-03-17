####################
# Auto_extend_nMP_range.R
# Extend getMetaPrograms to nMP 21:35, reusing existing geneNMF_outs.rds
# Skips any nMP value that already has a .rds in Metaprogrammes_Results/
# Run with gnmf env (needs GeneNMF package)
####################

library(GeneNMF)
library(RColorBrewer)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

# Load existing multiNMF output (already computed)
geneNMF.programs <- readRDS("geneNMF_outs.rds")

results_dir <- "Metaprogrammes_Results"
dir.create(results_dir, showWarnings = FALSE)

k_vals <- 21:35

for (k in k_vals) {
  rds_path <- file.path(results_dir, paste0("geneNMF_metaprograms_nMP_", k, ".rds"))

  # Skip if already generated
  if (file.exists(rds_path)) {
    message(paste0("nMP = ", k, " already exists, skipping"))
    next
  }

  message(paste0("Running getMetaPrograms with nMP = ", k))
  geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                          metric = "cosine",
                                          specificity.weight = 5,
                                          weight.explained = 0.5,
                                          nMP = k,
                                          min.confidence = 0.5)
  saveRDS(geneNMF.metaprograms, file = rds_path)

  # Save heatmap for each nMP
  n_colors <- min(k, 12)
  anno_colors <- brewer.pal(n = max(n_colors, 3), name = "Paired")
  anno_colors <- anno_colors[seq_len(length(geneNMF.metaprograms$metaprograms.genes))]
  names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

  png(file.path(results_dir, paste0("metaprograms_heatmap_nMP_", k, ".png")),
      width = 3000, height = 2500, res = 300)
  plotMetaPrograms(geneNMF.metaprograms,
                   annotation_colors = anno_colors,
                   similarity.cutoff = c(0, 1))
  dev.off()
  message(paste0("Saved nMP = ", k))
}

message("Extension complete. Run Auto_find_optimal_nmf.R to update metrics.")
