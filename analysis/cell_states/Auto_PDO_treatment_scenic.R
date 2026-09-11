####################
# Analysis registry:
#   Status: active terminal matched-treatment regulon workflow
#   Script: analysis/cell_states/Auto_PDO_treatment_scenic.R
#   Methodology: analysis/methodology/cell_states/Auto_PDO_treatment_scenic_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Runs SCENIC on four matched untreated/FLOT-treated PDO pairs and compares
#     regulon activity by treatment, sample, and paired patient direction.
#   Inputs:
#     - PDOs_outs/PDOs_merged.rds
#     - live hg38 cisTarget motif-ranking databases supplied by db_dir
#   Outputs:
#     - PDOs_outs/treatment_scenic/intermediate/*.rds
#     - PDOs_outs/treatment_scenic/tables/*.csv
#     - PDOs_outs/treatment_scenic/figures/*.pdf
#     - ephemeral PDOs_outs/treatment_scenic/cache/ SCENIC work files
#   Downstream use: none; persistent intermediates support replotting.
#   Cache/replot behavior: reuses SCENIC caches unless force=true; prepare_only
#     validates samples and writes the cell-count table without running SCENIC.
#   Run command: qsub run_treatment_scenic.sh
#   Conda env: dmtcp
####################

####################
# Auto_PDO_treatment_scenic.R
#
# Treatment-focused SCENIC workflow for matched FLOT pre/post PDO pairs.
# Runs SCENIC on the 8 matched samples (4 patients × Treated/Untreated)
# and compares regulon networks between Treated vs Untreated conditions,
# including patient-paired differential analysis.
#
# Input:
#   PDOs_outs/PDOs_merged.rds
#
# Output (in PDOs_outs/treatment_scenic/):
#   figures/  — regulon heatmaps, networks, volcano, paired-change plots
#   tables/   — differential regulon tables, network edges, per-patient stats
#   intermediate/ — SCENIC int/ folder, cached AUC/RSS matrices
#
# Env: dmtcp
# Usage:
#   Rscript analysis/cell_states/Auto_PDO_treatment_scenic.R
#   Rscript analysis/cell_states/Auto_PDO_treatment_scenic.R n_cores=12
#   Rscript analysis/cell_states/Auto_PDO_treatment_scenic.R prepare_only=true
####################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
library(data.table)
library(scales)
library(igraph)
library(ggraph)
library(tidygraph)
library(grid)
library(scales)
library(ggrepel)

## ------------------------------------------------------------------
## Arg parsing & helpers
## ------------------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x)) || !nzchar(x[1])) return(y)
  x[1]
}

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) next
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    out[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

to_flag <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) return(default)
  tolower(x) %in% c("true", "1", "yes", "y")
}

format_regulon_name <- function(x) {
  x <- gsub("_extended$", "", x)
  x <- gsub(" \\([0-9]+g\\)$", "", x)
  gsub(" \\([0-9]+ genes\\)$", "", x)
}

get_assay_matrix <- function(seurat_obj, slot_name = c("counts", "data")) {
  slot_name <- match.arg(slot_name)
  mat <- tryCatch(GetAssayData(seurat_obj, assay = "RNA", slot = slot_name), error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  mat <- tryCatch(LayerData(seurat_obj, assay = "RNA", layer = slot_name), error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  assay_obj <- seurat_obj@assays$RNA
  tryCatch(slot(assay_obj, slot_name),
           error = function(e) stop("Unable to retrieve RNA ", slot_name, " matrix from Seurat object."))
}

extract_regulon_targets <- function(x) {
  if (requireNamespace("GSEABase", quietly = TRUE) && methods::is(x, "GeneSet")) return(unique(GSEABase::geneIds(x)))
  if (is.character(x)) return(unique(x))
  if (is.list(x) && !is.null(x$gene)) return(unique(as.character(x$gene)))
  if (!is.null(names(x))) return(unique(names(x)))
  unique(as.character(x))
}

detect_db_files <- function(db_dir) {
  if (!dir.exists(db_dir)) stop("SCENIC database directory not found: ", db_dir)
  db_files <- list.files(db_dir, pattern = "\\.feather$", full.names = FALSE)
  db_files <- db_files[grepl("hg38|refseq-r80|hgnc", db_files, ignore.case = TRUE)]
  preferred <- db_files[grepl("mc9nr|refseq-r80", db_files, ignore.case = TRUE)]
  if (length(preferred) > 0) db_files <- preferred
  if (length(db_files) == 0) stop("No human cisTarget feather databases found in ", db_dir)
  unique(c(db_files[grepl("500bp", db_files, ignore.case = TRUE)][1],
           db_files[grepl("10kb", db_files, ignore.case = TRUE)][1],
           db_files))
}

patch_scenic_annotation_lookup <- function() {
  scenic_ns <- asNamespace("SCENIC")
  original_fun <- get("getDbAnnotations", envir = scenic_ns)
  patched_fun <- original_fun
  body(patched_fun) <- quote({
    dbAnnotFiles <- scenicOptions@settings$db_annotFiles
    if (!is.null(dbAnnotFiles)) {
      motifAnnotations <- NULL
      for (annotPath in dbAnnotFiles) {
        motifAnnot <- data.table::fread(annotPath)
        motifAnnot$annotationSource <- factor(motifAnnot$annotationSource)
        colnames(motifAnnot)[1] <- "motif"
        levels(motifAnnot$annotationSource) <- c(
          levels(motifAnnot$annotationSource),
          c("directAnnotation", "inferredBy_Orthology",
            "inferredBy_MotifSimilarity", "inferredBy_MotifSimilarity_n_Orthology"))
        motifAnnotations <- rbind(motifAnnotations, motifAnnot)
      }
    } else {
      if (is.na(getDatasetInfo(scenicOptions, "org"))) stop("Please provide an organism.")
      org <- getDatasetInfo(scenicOptions, "org")
      if (org == "hgnc") motifAnnotName <- "motifAnnotations_hgnc"
      if (org == "mgi") motifAnnotName <- "motifAnnotations_mgi"
      if (org == "dmel") motifAnnotName <- "motifAnnotations_dmel"
      if (!is.null(scenicOptions@settings$db_mcVersion) && scenicOptions@settings$db_mcVersion == "v8")
        motifAnnotName <- paste0(motifAnnotName, "_v8")
      annot_env <- new.env(parent = baseenv())
      data(list = motifAnnotName, package = "RcisTarget", envir = annot_env, verbose = FALSE)
      if (!exists(motifAnnotName, envir = annot_env, inherits = FALSE)) {
        v9_name <- paste0(motifAnnotName, "_v9")
        data(list = v9_name, package = "RcisTarget", envir = annot_env, verbose = FALSE)
        if (exists(v9_name, envir = annot_env, inherits = FALSE))
          assign(motifAnnotName, get(v9_name, envir = annot_env), envir = annot_env)
      }
      motifAnnotations <- get(motifAnnotName, envir = annot_env, inherits = FALSE)
    }
    return(motifAnnotations)
  })
  unlockBinding("getDbAnnotations", scenic_ns)
  assign("getDbAnnotations", patched_fun, envir = scenic_ns)
  lockBinding("getDbAnnotations", scenic_ns)
  invisible(TRUE)
}

patch_scenic_gene_filtering <- function() {
  scenic_ns <- asNamespace("SCENIC")
  original_fun <- get("geneFiltering", envir = scenic_ns)
  patched_fun <- original_fun
  body(patched_fun) <- quote({
    outFile_genesKept <- NULL
    dbFilePath <- NULL
    if (class(scenicOptions) == "ScenicOptions") {
      dbFilePath <- getDatabases(scenicOptions)[[1]]
      outFile_genesKept <- getIntName(scenicOptions, "genesKept")
    } else {
      dbFilePath <- scenicOptions[["dbFilePath"]]
      outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
    }
    if (is.null(dbFilePath)) stop("dbFilePath")
    if (is.data.frame(exprMat)) stop("data.frame expression matrices are not supported")
    if (any(table(rownames(exprMat)) > 1)) stop("Expression matrix rownames should be unique")
    if (inherits(exprMat, "Matrix") || inherits(exprMat, "sparseMatrix")) {
      nCountsPerGene <- Matrix::rowSums(exprMat, na.rm = TRUE)
      nCellsPerGene <- Matrix::rowSums(exprMat > 0, na.rm = TRUE)
    } else {
      nCountsPerGene <- rowSums(exprMat, na.rm = TRUE)
      nCellsPerGene <- rowSums(exprMat > 0, na.rm = TRUE)
    }
    genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
    nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
    genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
    library(RcisTarget)
    motifRankings <- importRankings(dbFilePath)
    genesInDatabase <- colnames(getRanking(motifRankings))
    genesKept <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
    if (!is.null(outFile_genesKept)) saveRDS(genesKept, file = outFile_genesKept)
    return(genesKept)
  })
  unlockBinding("geneFiltering", scenic_ns)
  assign("geneFiltering", patched_fun, envir = scenic_ns)
  lockBinding("geneFiltering", scenic_ns)
  invisible(TRUE)
}

scenic_gene_filtering_sparse <- function(exprMat, scenicOptions, minCountsPerGene, minSamples) {
  outFile_genesKept <- NULL
  dbFilePath <- NULL
  if (class(scenicOptions) == "ScenicOptions") {
    dbFilePath <- getDatabases(scenicOptions)[[1]]
    outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  } else {
    dbFilePath <- scenicOptions[["dbFilePath"]]
    outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
  }
  if (is.null(dbFilePath)) stop("dbFilePath")
  if (is.data.frame(exprMat)) stop("data.frame expression matrices are not supported")
  if (any(table(rownames(exprMat)) > 1)) stop("Expression matrix rownames should be unique")
  if (inherits(exprMat, "Matrix") || inherits(exprMat, "sparseMatrix")) {
    nCountsPerGene <- Matrix::rowSums(exprMat, na.rm = TRUE)
    nCellsPerGene <- Matrix::rowSums(exprMat > 0, na.rm = TRUE)
  } else {
    nCountsPerGene <- rowSums(exprMat, na.rm = TRUE)
    nCellsPerGene <- rowSums(exprMat > 0, na.rm = TRUE)
  }
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
  library(RcisTarget)
  motifRankings <- importRankings(dbFilePath)
  genesInDatabase <- colnames(getRanking(motifRankings))
  genesKept <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
  if (!is.null(outFile_genesKept)) saveRDS(genesKept, file = outFile_genesKept)
  genesKept
}

## ------------------------------------------------------------------
## Parse CLI arguments
## ------------------------------------------------------------------
arg_list <- parse_args(commandArgs(trailingOnly = TRUE))
prepare_only <- to_flag(arg_list[["prepare_only"]], default = FALSE)
n_cores <- as.integer(arg_list[["n_cores"]] %||% "8")
top_regulons_heatmap <- as.integer(arg_list[["top_regulons"]] %||% "8")
fdr_threshold <- as.numeric(arg_list[["fdr_threshold"]] %||% "0.05")
db_dir <- arg_list[["db_dir"]] %||% Sys.getenv("SCENIC_DB_DIR", unset = "")
if (!nzchar(db_dir)) db_dir <- "/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/cistarget_databases_rcistarget_mc9nr"
db_dir <- normalizePath(db_dir, winslash = "/", mustWork = FALSE)

## ------------------------------------------------------------------
## Matched sample definitions
## ------------------------------------------------------------------
patient_order <- c("SUR1070", "SUR1072", "SUR1090", "SUR1181")
matched_samples <- as.vector(rbind(
  paste0(patient_order, "_Untreated_PDO"),
  paste0(patient_order, "_Treated_PDO")
))
patient_cols <- c(SUR1070 = "#4C78A8", SUR1072 = "#59A14F",
                  SUR1090 = "#B07AA1", SUR1181 = "#F28E2B")
treatment_cols <- c(Untreated = "#D58B2D", Treated = "#374151")

## ------------------------------------------------------------------
## Directory structure
## ------------------------------------------------------------------
live_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
ephemeral_root <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"

out_dir <- file.path(live_root, "treatment_scenic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
tab_dir <- file.path(out_dir, "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

int_dir <- file.path(ephemeral_root, "treatment_scenic")
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(int_dir, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

message("=== Auto_PDO_treatment_scenic.R ===")
message("Live outputs:        ", out_dir)
message("Ephemeral SCENIC WD: ", int_dir)

## ------------------------------------------------------------------
## Load and subset to matched samples
## ------------------------------------------------------------------
message("Loading PDOs_merged.rds ...")
pdos <- readRDS(file.path(live_root, "PDOs_merged.rds"))

available_samples <- unique(pdos$orig.ident)
missing_samples <- setdiff(matched_samples, available_samples)
if (length(missing_samples) > 0) {
  stop("Missing matched samples in PDOs_merged.rds: ", paste(missing_samples, collapse = ", "))
}

pdos <- subset(pdos, subset = orig.ident %in% matched_samples)
n_cells <- ncol(pdos)
message("Matched-sample subset: ", n_cells, " cells across ", length(matched_samples), " samples")

## Derive treatment and patient metadata
pdos$treatment <- ifelse(grepl("_Treated_", pdos$orig.ident), "Treated", "Untreated")
pdos$patient <- sub("_(Treated|Untreated)_PDO$", "", pdos$orig.ident)
pdos$treatment <- factor(pdos$treatment, levels = c("Untreated", "Treated"))
pdos$patient <- factor(pdos$patient, levels = patient_order)

## Cell summary
cell_summary <- pdos@meta.data %>%
  group_by(orig.ident, patient, treatment) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(patient, treatment)
message("Per-sample cell counts:")
print(as.data.frame(cell_summary))

write.csv(cell_summary, file.path(tab_dir, "Auto_treatment_scenic_cell_summary.csv"), row.names = FALSE)

if (prepare_only) {
  message("prepare_only mode — cell summary saved, exiting.")
  quit(save = "no")
}

## ------------------------------------------------------------------
## SCENIC dependency checks
## ------------------------------------------------------------------
required_pkgs <- c("SCENIC", "AUCell", "RcisTarget", "GENIE3", "doRNG", "doMC")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required SCENIC packages: ", paste(missing_pkgs, collapse = ", "))
}
suppressPackageStartupMessages({
  library(SCENIC)
  library(AUCell)
  library(RcisTarget)
  library(GENIE3)
})

patch_scenic_annotation_lookup()
patch_scenic_gene_filtering()
db_files <- detect_db_files(db_dir)

## ------------------------------------------------------------------
## Expression matrix for SCENIC
## ------------------------------------------------------------------
counts_mat <- get_assay_matrix(pdos, "counts")
if (!inherits(counts_mat, "dgCMatrix")) counts_mat <- as(counts_mat, "dgCMatrix")
message("Expression matrix: ", nrow(counts_mat), " genes × ", ncol(counts_mat), " cells")

## Change to ephemeral SCENIC working directory (int/ folder lives here)
old_wd <- getwd()
setwd(int_dir)
on.exit(setwd(old_wd), add = TRUE)

scenicOptions <- initializeScenic(
  org = "hgnc", dbDir = db_dir, dbs = db_files,
  datasetTitle = "Auto_PDO_treatment_scenic", nCores = n_cores
)

min_counts_per_gene <- max(3 * 0.01 * ncol(counts_mat), 20)
min_samples <- max(0.01 * ncol(counts_mat), 20)

genes_kept_path <- file.path("int", "1.1_genesKept.Rds")
if (file.exists(genes_kept_path)) {
  message("Reusing existing genesKept: ", genes_kept_path)
  genes_kept <- readRDS(genes_kept_path)
} else {
  genes_kept <- scenic_gene_filtering_sparse(
    counts_mat, scenicOptions = scenicOptions,
    minCountsPerGene = min_counts_per_gene, minSamples = min_samples
  )
}

expr_mat_filtered <- counts_mat[genes_kept, , drop = FALSE]

db_tfs <- tryCatch(getDbTfs(scenicOptions), error = function(e) character(0))
focus_genes <- unique(c(
  intersect(db_tfs, rownames(expr_mat_filtered)),
  rownames(expr_mat_filtered)
))

if (length(focus_genes) >= 500) {
  expr_mat_use <- expr_mat_filtered[focus_genes, , drop = FALSE]
} else {
  expr_mat_use <- expr_mat_filtered
}
if (!is.matrix(expr_mat_use)) expr_mat_use <- as.matrix(expr_mat_use)

message("SCENIC input: ", nrow(expr_mat_use), " genes × ", ncol(expr_mat_use), " cells")

## ------------------------------------------------------------------
## SCENIC pipeline (with live checkpoint reuse)
## ------------------------------------------------------------------
live_regulon_auc_path <- file.path(out_dir, "intermediate", "regulon_auc.rds")
live_regulons_path <- file.path(out_dir, "intermediate", "regulons.rds")

if (file.exists(live_regulon_auc_path) && file.exists(live_regulons_path)) {
  message("Loading SCENIC results directly from live intermediate cache...")
  regulon_auc <- readRDS(live_regulon_auc_path)
  regulons <- readRDS(live_regulons_path)
  auc_mat <- getAUC(regulon_auc)
  message("SCENIC loaded: ", nrow(auc_mat), " regulons scored across ", ncol(auc_mat), " cells")
} else {
  if (!file.exists(file.path("int", "1.2_corrMat.Rds"))) {
    message("Running correlation ...")
    runCorrelation(expr_mat_use, scenicOptions)
  } else {
    message("Reusing existing correlation matrix.")
  }

  if (!file.exists(file.path("int", "1.4_GENIE3_linkList.Rds"))) {
    message("Running GENIE3 ...")
    runGenie3(expr_mat_use, scenicOptions, resumePreviousRun = TRUE, nParts = 10)
  } else {
    message("Reusing existing GENIE3 network.")
  }

  if (!file.exists(file.path("int", "1.6_tfModules_asDF.Rds"))) {
    message("Building co-expression modules ...")
    scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  } else {
    message("Reusing existing TF modules.")
  }

  if (!(file.exists(file.path("int", "2.6_regulons_asGeneSet.Rds")) &&
        file.exists(file.path("int", "2.6_regulons_asIncidMat.Rds")))) {
    message("Creating regulons ...")
    scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
  } else {
    message("Reusing existing regulons.")
  }

  if (!file.exists(file.path("int", "3.4_regulonAUC.Rds"))) {
    message("Scoring cells with AUCell ...")
    scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat = counts_mat)
  } else {
    message("Reusing existing regulon AUC scores.")
  }

  ## ------------------------------------------------------------------
  ## Load SCENIC results
  ## ------------------------------------------------------------------
  regulon_auc <- loadInt(scenicOptions, "aucell_regulonAUC")
  regulons <- loadInt(scenicOptions, "regulons")
  auc_mat <- getAUC(regulon_auc)
  message("SCENIC complete: ", nrow(auc_mat), " regulons scored across ", ncol(auc_mat), " cells")

  ## Cache core results to live
  saveRDS(regulon_auc, live_regulon_auc_path |>
            (\(x) { dir.create(dirname(x), recursive = TRUE, showWarnings = FALSE); x })())
  saveRDS(regulons, live_regulons_path |>
            (\(x) { dir.create(dirname(x), recursive = TRUE, showWarnings = FALSE); x })())
}

## ------------------------------------------------------------------
## Treatment-level mean AUC & RSS
## ------------------------------------------------------------------
treatment_map <- setNames(as.character(pdos$treatment), Cells(pdos))
treatment_map <- treatment_map[colnames(auc_mat)]
treatment_levels <- c("Untreated", "Treated")

mean_auc_treatment <- sapply(treatment_levels, function(trt) {
  cells <- names(treatment_map)[treatment_map == trt]
  rowMeans(auc_mat[, cells, drop = FALSE], na.rm = TRUE)
})
mean_auc_treatment <- as.matrix(mean_auc_treatment)

rss_treatment <- tryCatch(
  calcRSS(AUC = auc_mat, cellAnnotation = treatment_map),
  error = function(e) NULL
)
if (is.null(rss_treatment)) rss_treatment <- mean_auc_treatment
rss_treatment <- as.matrix(rss_treatment)
rss_treatment <- rss_treatment[, treatment_levels, drop = FALSE]

saveRDS(rss_treatment, file.path(out_dir, "intermediate", "rss_treatment.rds"))

## ------------------------------------------------------------------
## Per-patient mean AUC (for paired analysis)
## ------------------------------------------------------------------
sample_map <- setNames(as.character(pdos$orig.ident), Cells(pdos))
sample_map <- sample_map[colnames(auc_mat)]

mean_auc_sample <- sapply(matched_samples, function(s) {
  cells <- names(sample_map)[sample_map == s]
  if (length(cells) == 0) return(rep(NA_real_, nrow(auc_mat)))
  rowMeans(auc_mat[, cells, drop = FALSE], na.rm = TRUE)
})
mean_auc_sample <- as.matrix(mean_auc_sample)

rss_sample <- tryCatch(
  calcRSS(AUC = auc_mat, cellAnnotation = sample_map),
  error = function(e) NULL
)
if (is.null(rss_sample)) rss_sample <- mean_auc_sample
rss_sample <- as.matrix(rss_sample)
rss_sample <- rss_sample[, matched_samples, drop = FALSE]

saveRDS(mean_auc_sample, file.path(out_dir, "intermediate", "mean_auc_sample.rds"))
saveRDS(rss_sample, file.path(out_dir, "intermediate", "rss_sample.rds"))

## ------------------------------------------------------------------
## Differential regulon analysis: Treated vs Untreated
## ------------------------------------------------------------------
message("Computing differential regulon activity ...")

## Per-cell Wilcoxon test
diff_results <- bind_rows(lapply(rownames(auc_mat), function(regulon_id) {
  treated_vals <- auc_mat[regulon_id, names(treatment_map)[treatment_map == "Treated"]]
  untreated_vals <- auc_mat[regulon_id, names(treatment_map)[treatment_map == "Untreated"]]
  wt <- tryCatch(
    wilcox.test(treated_vals, untreated_vals, alternative = "two.sided", exact = FALSE),
    error = function(e) list(statistic = NA_real_, p.value = NA_real_)
  )
  data.frame(
    regulon = regulon_id,
    regulon_label = format_regulon_name(regulon_id),
    mean_auc_treated = mean(treated_vals, na.rm = TRUE),
    mean_auc_untreated = mean(untreated_vals, na.rm = TRUE),
    log2fc = log2((mean(treated_vals, na.rm = TRUE) + 1e-8) / (mean(untreated_vals, na.rm = TRUE) + 1e-8)),
    delta_auc = mean(treated_vals, na.rm = TRUE) - mean(untreated_vals, na.rm = TRUE),
    rss_treated = rss_treatment[regulon_id, "Treated"],
    rss_untreated = rss_treatment[regulon_id, "Untreated"],
    wilcox_p = wt$p.value,
    n_targets = length(extract_regulon_targets(regulons[[regulon_id]])),
    stringsAsFactors = FALSE
  )
}))

diff_results$wilcox_fdr <- p.adjust(diff_results$wilcox_p, method = "BH")
diff_results$direction <- ifelse(diff_results$delta_auc > 0, "up_in_treated",
                                  ifelse(diff_results$delta_auc < 0, "up_in_untreated", "no_change"))
diff_results$significant <- diff_results$wilcox_fdr < fdr_threshold
diff_results <- diff_results %>% arrange(wilcox_fdr, desc(abs(delta_auc)))

write.csv(diff_results, file.path(tab_dir, "Auto_treatment_scenic_diff_regulons.csv"), row.names = FALSE)

n_sig_up <- sum(diff_results$significant & diff_results$direction == "up_in_treated", na.rm = TRUE)
n_sig_dn <- sum(diff_results$significant & diff_results$direction == "up_in_untreated", na.rm = TRUE)
message("Significant regulons (FDR < ", fdr_threshold, "): ",
        n_sig_up, " up in treated, ", n_sig_dn, " up in untreated")

## ------------------------------------------------------------------
## Per-patient paired analysis
## ------------------------------------------------------------------
message("Computing per-patient paired regulon changes ...")

paired_results <- bind_rows(lapply(patient_order, function(patient) {
  untreated_sample <- paste0(patient, "_Untreated_PDO")
  treated_sample <- paste0(patient, "_Treated_PDO")
  bind_rows(lapply(rownames(auc_mat), function(regulon_id) {
    data.frame(
      patient = patient,
      regulon = regulon_id,
      regulon_label = format_regulon_name(regulon_id),
      mean_auc_untreated = mean_auc_sample[regulon_id, untreated_sample],
      mean_auc_treated = mean_auc_sample[regulon_id, treated_sample],
      delta_auc = mean_auc_sample[regulon_id, treated_sample] - mean_auc_sample[regulon_id, untreated_sample],
      stringsAsFactors = FALSE
    )
  }))
}))

## Paired sign-consistency: how many of the 4 patients show the same direction?
paired_consistency <- paired_results %>%
  group_by(regulon, regulon_label) %>%
  summarise(
    n_up = sum(delta_auc > 0, na.rm = TRUE),
    n_down = sum(delta_auc < 0, na.rm = TRUE),
    mean_delta = mean(delta_auc, na.rm = TRUE),
    sd_delta = sd(delta_auc, na.rm = TRUE),
    consistent_direction = ifelse(n_up >= 3, "up_in_treated",
                                   ifelse(n_down >= 3, "up_in_untreated", "inconsistent")),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(mean_delta)))

## Merge per-cell and paired results
diff_results_full <- diff_results %>%
  left_join(
    paired_consistency %>%
      select(regulon, n_up_patients = n_up, n_down_patients = n_down,
             paired_mean_delta = mean_delta, paired_sd_delta = sd_delta,
             paired_direction = consistent_direction),
    by = "regulon"
  )

write.csv(diff_results_full, file.path(tab_dir, "Auto_treatment_scenic_diff_regulons_with_paired.csv"), row.names = FALSE)
write.csv(paired_results, file.path(tab_dir, "Auto_treatment_scenic_paired_per_patient.csv"), row.names = FALSE)

## ------------------------------------------------------------------
## PLOT 1: Volcano plot (delta AUC vs -log10 FDR)
## ------------------------------------------------------------------
message("Generating plots ...")

volcano_df <- diff_results_full %>%
  mutate(
    neg_log10_fdr = -log10(pmax(wilcox_fdr, 1e-300)),
    color = case_when(
      significant & direction == "up_in_treated" ~ "Up in Treated",
      significant & direction == "up_in_untreated" ~ "Up in Untreated",
      TRUE ~ "Not significant"
    ),
    label = ifelse(significant & (rank(wilcox_fdr) <= 20 | abs(delta_auc) > quantile(abs(delta_auc), 0.95, na.rm = TRUE)),
                   regulon_label, NA_character_)
  )

volcano_cols <- c("Up in Treated" = "#374151", "Up in Untreated" = "#D58B2D", "Not significant" = "grey75")

p_volcano <- ggplot(volcano_df, aes(x = delta_auc, y = neg_log10_fdr, color = color)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf, na.rm = TRUE, color = "black") +
  scale_color_manual(values = volcano_cols) +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  labs(
    title = "Differential regulon activity: Treated vs Untreated",
    subtitle = paste0("Wilcoxon test; FDR < ", fdr_threshold,
                      " | ", n_sig_up, " up in treated, ", n_sig_dn, " up in untreated"),
    x = "Delta mean AUC (Treated - Untreated)",
    y = expression(-log[10](FDR)),
    color = NULL
  ) +
  theme_classic(base_size = 12)

ggsave(file.path(fig_dir, "Auto_treatment_scenic_volcano.pdf"), p_volcano, width = 12, height = 8)

## ------------------------------------------------------------------
## PLOT 2: Per-sample RSS heatmap (top regulons per treatment)
## ------------------------------------------------------------------
top_per_treatment <- unique(unlist(lapply(treatment_levels, function(trt) {
  vals <- sort(rss_treatment[, trt], decreasing = TRUE)
  names(vals)[seq_len(min(top_regulons_heatmap, length(vals)))]
})))
top_per_treatment <- top_per_treatment[!is.na(top_per_treatment)]

# Order matched samples: untreated first, then treated
ordered_samples <- c(
  paste0(patient_order, "_Untreated_PDO"),
  paste0(patient_order, "_Treated_PDO")
)

plot_rss_mat <- rss_sample[top_per_treatment, ordered_samples, drop = FALSE]
plot_rss_scaled <- t(scale(t(plot_rss_mat)))
plot_rss_scaled[!is.finite(plot_rss_scaled)] <- 0
rownames(plot_rss_scaled) <- format_regulon_name(rownames(plot_rss_scaled))

rss_col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

sample_treatment <- ifelse(grepl("_Treated_", ordered_samples), "Treated", "Untreated")
sample_patient <- sub("_(Treated|Untreated)_PDO$", "", ordered_samples)

ha_sample_rss <- HeatmapAnnotation(
  Treatment = factor(sample_treatment, levels = c("Untreated", "Treated")),
  Patient = factor(sample_patient, levels = patient_order),
  col = list(
    Treatment = treatment_cols,
    Patient = patient_cols
  ),
  show_annotation_name = TRUE
)

pdf(file.path(fig_dir, "Auto_treatment_scenic_rss_heatmap.pdf"), width = 10, height = 12, useDingbats = FALSE)
draw(
  Heatmap(
    plot_rss_scaled,
    name = "Scaled\nRSS",
    col = rss_col_fun,
    top_annotation = ha_sample_rss,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 12),
    column_names_rot = 90,
    heatmap_legend_param = list(title = "Scaled\nRSS")
  ),
  merge_legend = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
grid.text(
  "SCENIC regulon specificity: Treated vs Untreated",
  x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"),
  just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

## ------------------------------------------------------------------
## PLOT 3: Per-sample heatmap (all 8 samples)
## ------------------------------------------------------------------
## Top regulons: significant + consistent across patients
sig_consistent <- diff_results_full %>%
  filter(significant | (!is.na(paired_direction) & paired_direction != "inconsistent")) %>%
  arrange(wilcox_fdr) %>%
  head(60) %>%
  pull(regulon)

## Guarantee minimum from each direction
top_up <- diff_results_full %>%
  filter(direction == "up_in_treated") %>%
  arrange(wilcox_fdr) %>%
  head(15) %>%
  pull(regulon)
top_dn <- diff_results_full %>%
  filter(direction == "up_in_untreated") %>%
  arrange(wilcox_fdr) %>%
  head(15) %>%
  pull(regulon)

heatmap_regulons <- unique(c(sig_consistent, top_up, top_dn))
heatmap_regulons <- heatmap_regulons[heatmap_regulons %in% rownames(mean_auc_sample)]

if (length(heatmap_regulons) > 0) {
  # Order matched samples: untreated first, then treated
  ordered_samples <- c(
    paste0(patient_order, "_Untreated_PDO"),
    paste0(patient_order, "_Treated_PDO")
  )
  
  sample_heatmap_mat <- mean_auc_sample[heatmap_regulons, ordered_samples, drop = FALSE]
  
  # Center within each patient pair (subtract the pair's mean) to preserve magnitude
  sample_heatmap_scaled <- sample_heatmap_mat
  for (p in patient_order) {
    cols <- c(paste0(p, "_Untreated_PDO"), paste0(p, "_Treated_PDO"))
    pair_means <- rowMeans(sample_heatmap_mat[, cols, drop = FALSE])
    sample_heatmap_scaled[, cols[1]] <- sample_heatmap_mat[, cols[1]] - pair_means
    sample_heatmap_scaled[, cols[2]] <- sample_heatmap_mat[, cols[2]] - pair_means
  }
  sample_heatmap_scaled[!is.finite(sample_heatmap_scaled)] <- 0
  rownames(sample_heatmap_scaled) <- format_regulon_name(rownames(sample_heatmap_scaled))

  sample_treatment <- ifelse(grepl("_Treated_", ordered_samples), "Treated", "Untreated")
  sample_patient <- sub("_(Treated|Untreated)_PDO$", "", ordered_samples)

  ha_sample <- HeatmapAnnotation(
    Treatment = factor(sample_treatment, levels = c("Untreated", "Treated")),
    Patient = factor(sample_patient, levels = patient_order),
    col = list(
      Treatment = treatment_cols,
      Patient = patient_cols
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "left"
  )
  
  max_abs <- max(abs(sample_heatmap_scaled), na.rm = TRUE)
  auc_col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("#2166AC", "white", "#B2182B"))

  pdf(file.path(fig_dir, "Auto_treatment_scenic_sample_heatmap.pdf"), width = 12, height = 14, useDingbats = FALSE)
  draw(
    Heatmap(
      sample_heatmap_scaled,
      name = "Patient-scaled\nAUC",
      col = auc_col_fun,
      top_annotation = ha_sample,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      column_split = factor(sample_treatment, levels = c("Untreated", "Treated")),
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 9),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Scaled\nmean AUC")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid.text(
    "Regulon activity across matched FLOT pairs",
    x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"),
    just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()
}

## ------------------------------------------------------------------
## PLOT 4: Paired patient change plot (top significant regulons)
## ------------------------------------------------------------------
top_paired_regulons <- diff_results_full %>%
  filter(significant) %>%
  arrange(wilcox_fdr) %>%
  head(30) %>%
  pull(regulon)

if (length(top_paired_regulons) > 0) {
  paired_plot_df <- paired_results %>%
    filter(regulon %in% top_paired_regulons) %>%
    pivot_longer(cols = c(mean_auc_untreated, mean_auc_treated),
                 names_to = "condition", values_to = "mean_auc") %>%
    mutate(
      treatment = ifelse(condition == "mean_auc_treated", "Treated", "Untreated"),
      treatment = factor(treatment, levels = c("Untreated", "Treated")),
      regulon_label = factor(regulon_label, levels = unique(
        diff_results_full$regulon_label[diff_results_full$regulon %in% top_paired_regulons]
      ))
    )

  p_paired <- ggplot(paired_plot_df, aes(x = treatment, y = mean_auc, group = patient)) +
    geom_line(aes(color = patient), linewidth = 0.7, alpha = 0.8) +
    geom_point(aes(color = patient, shape = treatment), size = 2) +
    facet_wrap(~ regulon_label, scales = "free_y", ncol = 6) +
    scale_color_manual(values = patient_cols) +
    scale_shape_manual(values = c(Untreated = 16, Treated = 17)) +
    labs(
      title = "Patient-paired regulon activity changes (top significant)",
      x = NULL, y = "Mean AUCell score",
      color = "Patient", shape = "Treatment"
    ) +
    theme_classic(base_size = 10) +
    theme(
      strip.text = element_text(size = 7, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(file.path(fig_dir, "Auto_treatment_scenic_paired_changes.pdf"),
         p_paired, width = 18, height = max(6, 3 * ceiling(length(top_paired_regulons) / 6)))
}

## ------------------------------------------------------------------
## PLOT 5: Treatment-level regulon network
## ------------------------------------------------------------------
## Build bipartite network: top RSS regulons connected to Treated/Untreated
n_top_network <- 50

top_treated <- names(sort(rss_treatment[, "Treated"], decreasing = TRUE))[1:min(n_top_network, nrow(rss_treatment))]
top_untreated <- names(sort(rss_treatment[, "Untreated"], decreasing = TRUE))[1:min(n_top_network, nrow(rss_treatment))]
network_regulons <- unique(c(top_treated, top_untreated))

network_edge_df <- bind_rows(lapply(treatment_levels, function(trt) {
  vals <- rss_treatment[network_regulons, trt]
  vals <- vals[is.finite(vals) & vals > 0]
  threshold <- quantile(vals, 0.5, na.rm = TRUE)
  keep <- names(vals)[vals >= threshold]
  if (length(keep) == 0) return(NULL)
  data.frame(
    regulon_label = format_regulon_name(keep),
    treatment = trt,
    weight = as.numeric(vals[keep]),
    stringsAsFactors = FALSE
  )
})) %>% distinct(regulon_label, treatment, .keep_all = TRUE)

## Flag regulons shared between conditions
reg_degree <- network_edge_df %>% count(regulon_label, name = "n_treatments")

network_node_df <- data.frame(
  name = unique(c(network_edge_df$regulon_label, network_edge_df$treatment)),
  stringsAsFactors = FALSE
) %>%
  mutate(
    node_type = ifelse(name %in% treatment_levels, "Treatment", "Regulon"),
    node_group = case_when(
      name == "Treated" ~ "Treated",
      name == "Untreated" ~ "Untreated",
      TRUE ~ "Regulon"
    )
  ) %>%
  left_join(reg_degree, by = c("name" = "regulon_label")) %>%
  mutate(is_shared = !is.na(n_treatments) & n_treatments > 1)

## Color regulons by their differential direction
regulon_direction <- setNames(diff_results_full$direction, diff_results_full$regulon_label)
network_node_df$direction <- ifelse(
  network_node_df$node_type == "Regulon",
  regulon_direction[network_node_df$name],
  NA_character_
)
network_node_df$node_color <- case_when(
  network_node_df$name == "Treated" ~ "Treated",
  network_node_df$name == "Untreated" ~ "Untreated",
  network_node_df$is_shared ~ "Shared",
  network_node_df$direction == "up_in_treated" ~ "Treated-specific",
  network_node_df$direction == "up_in_untreated" ~ "Untreated-specific",
  TRUE ~ "Other"
)

network_graph <- tbl_graph(
  nodes = network_node_df,
  edges = network_edge_df %>% transmute(from = regulon_label, to = treatment, weight = weight),
  directed = FALSE
)

network_fill <- c(
  "Treated" = "#374151", "Untreated" = "#D58B2D",
  "Treated-specific" = "#6B7280", "Untreated-specific" = "#F0C27A",
  "Shared" = "#8B5CF6", "Other" = "grey70"
)

pdf(file.path(fig_dir, "Auto_treatment_scenic_network.pdf"), width = 18, height = 12, useDingbats = FALSE)
print(
  ggraph(network_graph, layout = "stress") +
    geom_edge_link(aes(width = weight, alpha = weight), colour = "grey70") +
    scale_edge_width(range = c(0.3, 2.0)) +
    scale_edge_alpha(range = c(0.2, 0.8)) +
    geom_node_point(
      aes(fill = node_color, shape = node_type,
          size = ifelse(node_type == "Treatment", 10, ifelse(is_shared, 5, 3.5))),
      colour = "black", stroke = 0.3
    ) +
    scale_size_identity() +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.5, max.overlaps = 40) +
    scale_shape_manual(values = c(Treatment = 21, Regulon = 22)) +
    scale_fill_manual(values = network_fill, drop = FALSE) +
    theme_void(base_size = 12) +
    labs(title = "SCENIC regulon network: Treated vs Untreated matched FLOT pairs") +
    guides(edge_width = "none", edge_alpha = "none", shape = "none",
           fill = guide_legend(title = "Node type"))
)
dev.off()

write.csv(network_edge_df, file.path(tab_dir, "Auto_treatment_scenic_network_edges.csv"), row.names = FALSE)

## ------------------------------------------------------------------
## PLOT 6: All-regulon clustering heatmap (MPs hierarchically clustered)
## ------------------------------------------------------------------
global_top_regulons <- names(sort(apply(mean_auc_sample[, matched_samples, drop = FALSE], 1, max), decreasing = TRUE))
global_top_regulons <- global_top_regulons[1:min(100, length(global_top_regulons))]
global_top_regulons <- global_top_regulons[!is.na(global_top_regulons)]

if (length(global_top_regulons) > 0) {
  clust_mat <- mean_auc_sample[global_top_regulons, matched_samples, drop = FALSE]
  clust_scaled <- t(scale(t(clust_mat)))
  clust_scaled[!is.finite(clust_scaled)] <- 0
  rownames(clust_scaled) <- format_regulon_name(rownames(clust_scaled))

  matched_sample_treatment <- ifelse(grepl("_Treated_", matched_samples), "Treated", "Untreated")
  matched_sample_patient <- sub("_(Treated|Untreated)_PDO$", "", matched_samples)

  ha_matched_sample <- HeatmapAnnotation(
    Treatment = factor(matched_sample_treatment, levels = c("Untreated", "Treated")),
    Patient = factor(matched_sample_patient, levels = patient_order),
    col = list(
      Treatment = treatment_cols,
      Patient = patient_cols
    ),
    show_annotation_name = TRUE,
    annotation_name_side = "left"
  )

  pdf(file.path(fig_dir, "Auto_treatment_scenic_clustering_heatmap.pdf"), width = 12, height = 14, useDingbats = FALSE)
  draw(
    Heatmap(
      clust_scaled,
      name = "Scaled\nmean AUC",
      col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
      top_annotation = ha_matched_sample,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 9),
      column_names_rot = 45,
      heatmap_legend_param = list(title = "Scaled\nmean AUC")
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid.text(
    "Top 100 regulons: unsupervised clustering across matched pairs",
    x = unit(4, "mm"), y = unit(1, "npc") - unit(4, "mm"),
    just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")
  )
  dev.off()
}

## ------------------------------------------------------------------
## PLOT 7: Regulon target summary for top differential regulons
## ------------------------------------------------------------------
top_diff <- diff_results_full %>%
  filter(significant) %>%
  arrange(wilcox_fdr) %>%
  head(50)

if (nrow(top_diff) > 0) {
  regulon_target_df <- bind_rows(lapply(seq_len(nrow(top_diff)), function(i) {
    reg_name <- top_diff$regulon[i]
    reg_targets <- extract_regulon_targets(regulons[[reg_name]])
    data.frame(
      regulon = reg_name,
      regulon_label = top_diff$regulon_label[i],
      direction = top_diff$direction[i],
      wilcox_fdr = top_diff$wilcox_fdr[i],
      delta_auc = top_diff$delta_auc[i],
      n_targets = length(reg_targets),
      n_up_patients = top_diff$n_up_patients[i],
      n_down_patients = top_diff$n_down_patients[i],
      targets_preview = paste(head(reg_targets, 50), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(regulon_target_df, file.path(tab_dir, "Auto_treatment_scenic_top_regulon_targets.csv"), row.names = FALSE)
}

## ------------------------------------------------------------------
## Summary
## ------------------------------------------------------------------
summary_df <- data.frame(
  n_matched_samples = length(matched_samples),
  n_patients = length(patient_order),
  n_cells_total = n_cells,
  n_regulons = nrow(auc_mat),
  n_sig_up_treated = n_sig_up,
  n_sig_up_untreated = n_sig_dn,
  fdr_threshold = fdr_threshold,
  n_cores = n_cores,
  db_dir = db_dir,
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(tab_dir, "Auto_treatment_scenic_summary.csv"), row.names = FALSE)

## ------------------------------------------------------------------
## State x Treatment analysis
## ------------------------------------------------------------------
message("Loading canonical centred states...")
final_states_path <- file.path(live_root, "centred_mp_refinement", "centred_refined_noreg_states.rds")
if (file.exists(final_states_path)) {
  final_states <- readRDS(final_states_path)
  
  # Ensure we only keep matched cells
  common_cells <- intersect(colnames(auc_mat), names(final_states))
  
  if (length(common_cells) > 0) {
    message("Found ", length(common_cells), " cells with canonical state assignments.")
    
    # Filter to discrete canonical states
    state_levels <- c("Classic proliferation", "Columnar-to-intestinal", 
                      "Glandular differentiation", "Stress-adaptive",
                      "ECM-remodelling", "Motile-cilia differentiation")
    
    canonical_cells <- common_cells[final_states[common_cells] %in% state_levels]
    
    auc_mat_canon <- auc_mat[, canonical_cells, drop = FALSE]
    canon_states <- factor(final_states[canonical_cells], levels = state_levels)
    canon_treats <- factor(treatment_map[canonical_cells], levels = c("Untreated", "Treated"))
    
    state_treat_map <- paste(canon_states, canon_treats, sep = "_")
    names(state_treat_map) <- canonical_cells
    
    # Reorder state_treat_map to have State(Untreated) next to State(Treated)
    st_levels <- paste(rep(state_levels, each = 2), rep(c("Untreated", "Treated"), length(state_levels)), sep = "_")
    state_treat_map <- factor(state_treat_map, levels = st_levels)
    
    message("Calculating mean AUC by State x Treatment...")
    mean_auc_st <- sapply(st_levels, function(g) {
      cells <- names(state_treat_map)[state_treat_map == g]
      if (length(cells) == 0) return(rep(0, nrow(auc_mat_canon)))
      rowMeans(auc_mat_canon[, cells, drop = FALSE], na.rm = TRUE)
    })
    rownames(mean_auc_st) <- rownames(auc_mat_canon)
    
    message("Finding top regulons per state for treated vs untreated differences...")
    # Find top differential regulons per state (Treated vs Untreated)
    top_diff_st_regulons <- unique(unlist(lapply(state_levels, function(s) {
      untreated_cells <- names(canon_states)[canon_states == s & canon_treats == "Untreated"]
      treated_cells <- names(canon_states)[canon_states == s & canon_treats == "Treated"]
      if(length(untreated_cells) > 5 && length(treated_cells) > 5) {
        pvals <- sapply(rownames(auc_mat_canon), function(r) {
          suppressWarnings(wilcox.test(auc_mat_canon[r, treated_cells], auc_mat_canon[r, untreated_cells])$p.value)
        })
        pvals[is.na(pvals)] <- 1
        log2fc <- sapply(rownames(auc_mat_canon), function(r) {
          # add pseudocount to avoid inf
          log2((mean(auc_mat_canon[r, treated_cells]) + 1e-6) / (mean(auc_mat_canon[r, untreated_cells]) + 1e-6))
        })
        sig <- pvals < 0.05
        up <- names(sort(log2fc[sig & log2fc > 0], decreasing = TRUE))[1:min(5, sum(sig & log2fc > 0))]
        down <- names(sort(log2fc[sig & log2fc < 0], decreasing = FALSE))[1:min(5, sum(sig & log2fc < 0))]
        return(c(up, down))
      }
      return(NULL)
    })))
    top_diff_st_regulons <- top_diff_st_regulons[!is.na(top_diff_st_regulons)]
    
    if (length(top_diff_st_regulons) > 0) {
      plot_mat_st <- mean_auc_st[top_diff_st_regulons, , drop = FALSE]
      
      # 1. Z-score scaling (across row)
      plot_z_st <- t(scale(t(plot_mat_st)))
      plot_z_st[!is.finite(plot_z_st)] <- 0
      rownames(plot_z_st) <- format_regulon_name(rownames(plot_z_st))
      
      # 2. Paired scaling (subtract pair mean)
      plot_pair_st <- plot_mat_st
      for (s in state_levels) {
        cols <- paste0(s, c("_Untreated", "_Treated"))
        pair_means <- rowMeans(plot_mat_st[, cols, drop = FALSE])
        plot_pair_st[, cols[1]] <- plot_mat_st[, cols[1]] - pair_means
        plot_pair_st[, cols[2]] <- plot_mat_st[, cols[2]] - pair_means
      }
      plot_pair_st[!is.finite(plot_pair_st)] <- 0
      rownames(plot_pair_st) <- format_regulon_name(rownames(plot_pair_st))
      
      # Use full names with newlines for the annotation to prevent overlap
      newline_states <- c(
        "Classic proliferation" = "Classic\nproliferation",
        "Columnar-to-intestinal" = "Columnar-to-\nintestinal",
        "Glandular differentiation" = "Glandular\ndifferentiation",
        "Stress-adaptive" = "Stress-\nadaptive",
        "ECM-remodelling" = "ECM-\nremodelling",
        "Motile-cilia differentiation" = "Motile-cilia\ndifferentiation"
      )
      
      ha_st <- HeatmapAnnotation(
        Treatment = rep(c("Untreated", "Treated"), length(state_levels)),
        State = rep(newline_states[state_levels], each = 2),
        col = list(
          Treatment = treatment_cols,
          State = setNames(c("#E41A1C", "#4DAF4A", "#FF7F00", "#984EA3", "#A65628", "#F781BF"), newline_states[state_levels])
        ),
        show_annotation_name = TRUE
      )
      
      max_abs_st <- max(abs(plot_pair_st), na.rm = TRUE)
      col_fun_z <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
      col_fun_pair <- colorRamp2(c(-max_abs_st, 0, max_abs_st), c("#2166AC", "white", "#B2182B"))
      
      pdf(file.path(fig_dir, "Auto_treatment_scenic_state_treatment_heatmap.pdf"), width = 12, height = max(8, length(top_diff_st_regulons)*0.2), useDingbats = FALSE)
      
      # Page 1: Z-score scaled, Clustered rows
      draw(
        Heatmap(
          plot_z_st,
          name = "Z-score\nAUC",
          col = col_fun_z,
          top_annotation = ha_st,
          cluster_rows = TRUE,
          cluster_columns = FALSE, 
          column_split = rep(newline_states[state_levels], each = 2),
          show_column_names = FALSE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "Z-score\nAUC"),
          column_title_rot = 0,
          column_title_gp = gpar(fontsize = 10)
        ),
        merge_legend = TRUE,
        column_title = "Across-row Z-score Scaled (Clustered)"
      )
      
      # Page 2: Paired scaled, Clustered rows
      draw(
        Heatmap(
          plot_pair_st,
          name = "Pair-scaled\nAUC",
          col = col_fun_pair,
          top_annotation = ha_st,
          cluster_rows = TRUE, 
          cluster_columns = FALSE, 
          column_split = rep(newline_states[state_levels], each = 2),
          show_column_names = FALSE,
          row_names_side = "left",
          row_names_gp = gpar(fontsize = 8),
          heatmap_legend_param = list(title = "Pair-scaled\nAUC"),
          column_title_rot = 0,
          column_title_gp = gpar(fontsize = 10)
        ),
        merge_legend = TRUE,
        column_title = "State Pair-scaled (Clustered)"
      )
      
      dev.off()
      message("Saved state x treatment heatmap.")
      
      # Save the mean AUC per state x treatment table
      write.csv(mean_auc_st, file.path(tab_dir, "Auto_treatment_scenic_state_treatment_meanAUC.csv"), row.names = TRUE)
    } else {
      message("Not enough cells/significant differences per state for heatmap.")
    }
  } else {
    message("No cell overlap with canonical states.")
  }
} else {
  message("States file not found: ", final_states_path)
}

message("=== Complete ===")
message("Live outputs:        ", out_dir)
message("Ephemeral SCENIC WD: ", int_dir)
message("Total regulons:      ", nrow(auc_mat))
message("Significant (FDR < ", fdr_threshold, "): ",
        n_sig_up, " up treated, ", n_sig_dn, " up untreated")
