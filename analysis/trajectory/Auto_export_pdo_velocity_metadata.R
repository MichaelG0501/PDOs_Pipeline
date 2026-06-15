####################
# Auto_export_pdo_velocity_metadata.R
#
# Export PDO per-cell metadata and per-sample barcode lists for RNA velocity.
# Uses finalized states for UMAP coloring and the pre-unresolved-relabel
# four-state vector for state-direction summaries.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
})

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_velocity_PDO"
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "barcodes"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

message("Loading PDO object and state vectors...")
pdos <- readRDS("PDOs_merged.rds")
state_four <- readRDS("Auto_PDO_states_noreg.rds")
state_final <- readRDS("Auto_PDO_final_states.rds")

state_four_names <- names(state_four)
state_four <- as.character(state_four)
names(state_four) <- state_four_names

state_final_names <- names(state_final)
state_final <- as.character(state_final)
names(state_final) <- state_final_names

cells <- Cells(pdos)
umap <- Embeddings(pdos, "umap")
umap <- umap[cells, 1:2, drop = FALSE]

meta <- pdos@meta.data[cells, , drop = FALSE]
meta$cell_id <- cells
meta$sample <- as.character(meta$orig.ident)
meta$batch_type <- ifelse(grepl("_(Untreated|Treated)_PDO$", meta$sample), "New_batch", "Cynthia_batch")
meta$treatment <- dplyr::case_when(
  grepl("_Untreated_PDO$", meta$sample) ~ "Untreated",
  grepl("_Treated_PDO$", meta$sample) ~ "Treated",
  TRUE ~ "Cynthia"
)
meta$raw_barcode <- mapply(
  function(cell_id, sample) {
    prefix <- paste0(sample, "_")
    if (startsWith(cell_id, prefix)) {
      substring(cell_id, nchar(prefix) + 1)
    } else {
      sub("^.*_", "", cell_id)
    }
  },
  meta$cell_id,
  meta$sample,
  USE.NAMES = FALSE
)
meta$state_four <- state_four[meta$cell_id]
meta$state_final <- state_final[meta$cell_id]
meta$umap_1 <- umap[, 1]
meta$umap_2 <- umap[, 2]

keep_cols <- c(
  "cell_id", "sample", "raw_barcode", "batch_type", "treatment",
  "state_four", "state_final", "umap_1", "umap_2",
  "nCount_RNA", "nFeature_RNA", "SUR", "Batch"
)
keep_cols <- keep_cols[keep_cols %in% colnames(meta)]
meta_out <- meta[, keep_cols, drop = FALSE]

fwrite(meta_out, file.path(out_dir, "tables", "Auto_pdo_velocity_cell_metadata.csv"))

sample_names <- sort(unique(meta_out$sample))
for (sample_name in sample_names) {
  barcodes <- sort(unique(meta_out$raw_barcode[meta_out$sample == sample_name]))
  writeLines(barcodes, file.path(out_dir, "barcodes", paste0(sample_name, "_qc_barcodes.tsv")))
}

new_bam_root <- "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/cellranger"
cynthia_fastq_root <- "/rds/general/project/spatialtranscriptomics/ephemeral/PDOs"
wd <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"

manifest <- data.frame(
  sample = sample_names,
  batch_type = ifelse(grepl("_(Untreated|Treated)_PDO$", sample_names), "New_batch", "Cynthia_batch"),
  treatment = dplyr::case_when(
    grepl("_Untreated_PDO$", sample_names) ~ "Untreated",
    grepl("_Treated_PDO$", sample_names) ~ "Treated",
    TRUE ~ "Cynthia"
  ),
  stringsAsFactors = FALSE
)
manifest$pool <- ifelse(
  manifest$treatment == "Untreated", "PDOs_Untreated",
  ifelse(manifest$treatment == "Treated", "PDOs_Treated", NA_character_)
)
manifest$fastq_dir <- ifelse(
  manifest$batch_type == "Cynthia_batch",
  file.path(cynthia_fastq_root, manifest$sample),
  NA_character_
)
manifest$cellranger_out <- ifelse(
  manifest$batch_type == "Cynthia_batch",
  file.path(wd, "PDOs_outs", "Auto_velocity_PDO", "cellranger", manifest$sample, "outs"),
  file.path(new_bam_root, manifest$pool, "outs")
)
manifest$bam <- file.path(manifest$cellranger_out, "possorted_genome_bam.bam")
manifest$barcodes_file <- file.path(wd, "PDOs_outs", out_dir, "barcodes", paste0(manifest$sample, "_qc_barcodes.tsv"))
manifest$n_cells <- as.integer(table(factor(meta_out$sample, levels = manifest$sample)))
manifest$has_bam <- file.exists(manifest$bam)

fwrite(manifest, file.path(out_dir, "tables", "Auto_pdo_velocity_sample_manifest.csv"))

summary_tbl <- meta_out %>%
  count(sample, batch_type, treatment, state_four, name = "cells") %>%
  arrange(batch_type, sample, state_four)
fwrite(summary_tbl, file.path(out_dir, "tables", "Auto_pdo_velocity_state_four_summary.csv"))

message("Wrote velocity metadata for ", nrow(meta_out), " cells and ", nrow(manifest), " samples.")
