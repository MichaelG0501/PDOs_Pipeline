library(GeneNMF)
library(RColorBrewer)
library(msigdbr)
library(fgsea)
library(UCell)
library(Seurat)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

####################
# Load input data
####################
if (!file.exists("PDOs_list_PDOs.rds")) stop("Input file PDOs_list_PDOs.rds not found!")
pdos.list <- readRDS("PDOs_list_PDOs.rds")
pdos.list$SUR843T3_PDO <- NULL

####################
# Step 1: multiNMF
####################
if (file.exists("geneNMF_outs.rds")) {
  message("Loading existing geneNMF_outs.rds")
  geneNMF.programs <- readRDS("geneNMF_outs.rds")
} else {
  message("Running multiNMF...")
  geneNMF.programs <- multiNMF(pdos.list, assay="RNA", k=4:9, min.exp = 0.05)
  saveRDS(geneNMF.programs, file="geneNMF_outs.rds")
}

####################
# Step 2: Loop over nMP values (4:20)
####################
dir.create("Metaprogrammes_Results", showWarnings = FALSE)
k_vals <- 4:20
for (k in k_vals) {
  rds_path <- file.path("Metaprogrammes_Results", paste0("geneNMF_metaprograms_nMP_", k, ".rds"))
  png_path <- file.path("Metaprogrammes_Results", paste0("metaprograms_heatmap_nMP_", k, ".png"))
  
  if (file.exists(rds_path) && file.exists(png_path)) {
    message(paste0("nMP = ", k, " already exists, skipping."))
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

  png(png_path, width = 3000, height = 2500, res = 300)
  plotMetaPrograms(geneNMF.metaprograms,
                   annotation_colors = anno_colors,
                   similarity.cutoff = c(0, 1))
  dev.off()
  message(paste0("Saved nMP = ", k))
}

####################
# Step 3: Optimal nMP=13
####################
default_nMP <- 13
rds_default <- file.path("Metaprogrammes_Results", paste0("geneNMF_metaprograms_nMP_", default_nMP, ".rds"))

if (file.exists("MP_outs_default.rds")) {
  message("Loading existing MP_outs_default.rds")
  geneNMF.metaprograms <- readRDS("MP_outs_default.rds")
} else {
  if (!file.exists(rds_default)) stop(paste0("Default nMP file not found: ", rds_default))
  message(paste0("Setting nMP = ", default_nMP, " as default."))
  geneNMF.metaprograms <- readRDS(rds_default)
  saveRDS(geneNMF.metaprograms, file = "MP_outs_default.rds")
  
  anno_colors <- brewer.pal(n = min(default_nMP, 12), name = "Paired")
  anno_colors <- anno_colors[seq_len(length(geneNMF.metaprograms$metaprograms.genes))]
  names(anno_colors) <- names(geneNMF.metaprograms$metaprograms.genes)

  png("metaprograms_heatmap.png", width = 3000, height = 2500, res = 300)
  plotMetaPrograms(geneNMF.metaprograms,
                   annotation_colors = anno_colors,
                   similarity.cutoff = c(0, 1))
  dev.off()
}

####################
# Step 4: GSEA and UCell Scoring
####################
if (!file.exists("PDOs_merged.rds")) stop("Input file PDOs_merged.rds not found!")

if (file.exists("GO_outs.rds")) {
  message("Loading existing GO_outs.rds")
  top_p <- readRDS("GO_outs.rds")
} else {
  message("Running GSEA...")
  pdos <- readRDS("PDOs_merged.rds")
  pdos <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")
  
  # Filter MPs by silhouette and sample coverage as per AGENTS.md
  bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
  bad_mp_names <- paste0("MP", bad_mps)
  
  coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
  names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
  low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
  
  exclude_mps <- unique(c(bad_mp_names, low_coverage_mps))
  if (length(exclude_mps) > 0) {
    message(paste0("Excluding MPs: ", paste(exclude_mps, collapse=", ")))
  }
  
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  mp.genes <- mp.genes[!names(mp.genes) %in% exclude_mps]
  
  top_p <- lapply(mp.genes, function(program) {
    runGSEA(program, universe=rownames(pdos), category = "C5", subcategory = "GO:BP")
  })
  saveRDS(top_p, "GO_outs.rds")
}

if (file.exists("PDOs_final.rds")) {
  message("Loading existing PDOs_final.rds")
  pdos <- readRDS("PDOs_final.rds")
} else {
  message("Running UCell scoring...")
  if (!exists("pdos")) {
    pdos <- readRDS("PDOs_merged.rds")
    pdos <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")
  }
  
  # Re-calculate mp.genes in case pdos was just loaded
  bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
  bad_mp_names <- paste0("MP", bad_mps)
  coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
  names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
  low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
  exclude_mps <- unique(c(bad_mp_names, low_coverage_mps))
  
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
  mp.genes <- mp.genes[!names(mp.genes) %in% exclude_mps]
  
  pdos <- AddModuleScore_UCell(pdos, features = mp.genes, ncores=4, name = "")
  
  # Save UCell scores separately for lightweight downstream use
  ucell_cols <- names(mp.genes)
  ucell_scores <- pdos@meta.data[, ucell_cols, drop = FALSE]
  saveRDS(ucell_scores, "UCell_scores_filtered.rds")
  message("Saved UCell_scores_filtered.rds")
  
  saveRDS(pdos, "PDOs_final.rds")
}

####################
# Step 5: Plotting
####################
message("Generating final plots...")

# Ensure mp.genes is correctly defined for plotting
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
exclude_mps <- unique(c(bad_mp_names, low_coverage_mps))

mp.genes <- geneNMF.metaprograms$metaprograms.genes
mp.genes <- mp.genes[!names(mp.genes) %in% exclude_mps]

png("vln_origident.png", width = 5000, height = 2500, res = 300)
print(VlnPlot(pdos, features = names(mp.genes), group.by = "orig.ident", pt.size = 0, ncol = 5))
dev.off()

if ("Clinical.response.at.OG.MDT..responder.non.responder" %in% colnames(pdos@meta.data)) {
  png("vln_clinical_response.png", width = 3500, height = 2500, res = 300)
  print(VlnPlot(pdos, features = names(mp.genes), group.by = "Clinical.response.at.OG.MDT..responder.non.responder", pt.size = 0, ncol = 5))
  dev.off()
}

if ("Batch" %in% colnames(pdos@meta.data)) {
  png("vln_batch.png", width = 3500, height = 2500, res = 300)
  print(VlnPlot(pdos, features = names(mp.genes), group.by = "Batch", pt.size = 0, ncol = 5))
  dev.off()
}

message("geneNMF.R completed successfully.")
