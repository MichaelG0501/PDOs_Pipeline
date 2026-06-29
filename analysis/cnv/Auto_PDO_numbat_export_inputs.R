####################
# Auto_PDO_numbat_export_inputs.R
#
# Export one raw-count matrix and metadata map per PDO sample for Numbat.
# Numbat is run on raw 10x cell barcodes; downstream PDO objects use
# sample-prefixed cell IDs, so this script records both identifiers.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

out_dir <- "Auto_PDO_numbat"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "by_samples"), recursive = TRUE, showWarnings = FALSE)

velocity_manifest_path <- file.path(
  out_root,
  "Auto_velocity_PDO",
  "tables",
  "Auto_pdo_velocity_sample_manifest.csv"
)
if (!file.exists(velocity_manifest_path)) {
  stop("Missing velocity manifest used for BAM/barcode paths: ", velocity_manifest_path)
}

excluded_samples <- c("SUR843T3_PDO")
velocity_manifest <- fread(velocity_manifest_path)
velocity_manifest <- velocity_manifest[!sample %in% excluded_samples]

get_counts <- function(obj) {
  suppressWarnings({
    tryCatch(
      GetAssayData(obj, assay = "RNA", layer = "counts"),
      error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
    )
  })
}

sample_to_cell_id <- function(sample_name, raw_barcode) {
  paste(sample_name, raw_barcode, sep = "_")
}

resolve_bam <- function(sample_name, batch_type, manifest_bam) {
  cynthia_bam <- file.path(
    out_root,
    "Auto_velocity_PDO",
    "cellranger",
    sample_name,
    "outs",
    "possorted_genome_bam.bam"
  )
  coord_bam <- file.path(
    out_root,
    "Auto_velocity_PDO",
    "coord",
    paste0(sample_name, ".qc.coord.bam")
  )

  ####################
  # Live Cell Ranger BAM fallback for regenerated Numbat inputs. Older Cynthia
  # PDOs are sample-specific; newer treated/untreated PDOs use pooled Cell
  # Ranger BAMs with sample-specific barcode files.
  live_cellranger_root <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Cellranger_outs"
  live_sample_bam <- file.path(
    live_cellranger_root,
    sample_name,
    "outs",
    "possorted_genome_bam.bam"
  )
  live_pool_name <- case_when(
    grepl("_Treated_PDO$", sample_name) ~ "PDOs_Treated",
    grepl("_Untreated_PDO$", sample_name) ~ "PDOs_Untreated",
    TRUE ~ NA_character_
  )
  live_pool_bam <- if (!is.na(live_pool_name)) {
    file.path(live_cellranger_root, live_pool_name, "outs", "possorted_genome_bam.bam")
  } else {
    NA_character_
  }
  ####################

  candidates <- if (identical(batch_type, "Cynthia_batch")) {
    c(cynthia_bam, coord_bam, live_sample_bam, manifest_bam)
  } else {
    c(manifest_bam, coord_bam, live_pool_bam, live_sample_bam)
  }
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  normalizePath(hit[1], mustWork = TRUE)
}

manifest_rows <- list()
for (i in seq_len(nrow(velocity_manifest))) {
  sample_name <- velocity_manifest$sample[i]
  sample_rds <- file.path(out_root, "by_samples", sample_name, paste0(sample_name, ".rds"))
  if (!file.exists(sample_rds)) {
    warning("Skipping sample with no post-QC Seurat RDS: ", sample_name)
    next
  }

  sample_out <- file.path(out_root, out_dir, "by_samples", sample_name)
  input_dir <- file.path(sample_out, "input")
  numbat_dir <- file.path(sample_out, "numbat")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(numbat_dir, recursive = TRUE, showWarnings = FALSE)

  message("Exporting Numbat inputs for ", sample_name)
  obj <- readRDS(sample_rds)
  counts <- get_counts(obj)
  raw_barcodes <- colnames(counts)
  cell_ids <- sample_to_cell_id(sample_name, raw_barcodes)

  count_rds <- file.path(input_dir, paste0("Auto_", sample_name, "_counts_raw_barcodes.rds"))
  metadata_csv <- file.path(input_dir, paste0("Auto_", sample_name, "_cell_map.csv"))
  saveRDS(counts, count_rds)

  meta <- obj@meta.data
  meta$raw_barcode <- raw_barcodes
  meta$cell_id <- cell_ids
  meta$sample <- sample_name
  meta_out <- meta[, unique(c(
    "cell_id", "sample", "raw_barcode", "orig.ident", "nCount_RNA",
    "nFeature_RNA", "percent.mt", "Batch", "SUR"
  )[c(
    "cell_id", "sample", "raw_barcode", "orig.ident", "nCount_RNA",
    "nFeature_RNA", "percent.mt", "Batch", "SUR"
  ) %in% colnames(meta)]), drop = FALSE]
  fwrite(meta_out, metadata_csv)

  barcodes_file <- velocity_manifest$barcodes_file[i]
  if (!file.exists(barcodes_file)) {
    warning("Missing barcode file for ", sample_name, ": ", barcodes_file)
  } else {
    listed <- readLines(barcodes_file, warn = FALSE)
    missing_in_counts <- setdiff(listed, raw_barcodes)
    if (length(missing_in_counts) > 0) {
      warning(
        sample_name,
        " barcode file has ",
        length(missing_in_counts),
        " barcodes absent from the sample RDS counts."
      )
    }
  }

  bam <- resolve_bam(sample_name, velocity_manifest$batch_type[i], velocity_manifest$bam[i])
  allele_file <- file.path(sample_out, paste0(sample_name, "_allele_counts.tsv.gz"))
  clone_post_file <- file.path(numbat_dir, "clone_post_2.tsv")
  joint_post_file <- file.path(numbat_dir, "joint_post_2.tsv")

  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    sample = sample_name,
    batch_type = velocity_manifest$batch_type[i],
    treatment = velocity_manifest$treatment[i],
    pool = velocity_manifest$pool[i],
    bam = bam,
    barcodes_file = normalizePath(barcodes_file, mustWork = FALSE),
    count_rds = normalizePath(count_rds, mustWork = TRUE),
    metadata_csv = normalizePath(metadata_csv, mustWork = TRUE),
    sample_out_dir = normalizePath(sample_out, mustWork = TRUE),
    numbat_dir = normalizePath(numbat_dir, mustWork = TRUE),
    allele_file = normalizePath(allele_file, mustWork = FALSE),
    clone_post_file = normalizePath(clone_post_file, mustWork = FALSE),
    joint_post_file = normalizePath(joint_post_file, mustWork = FALSE),
    n_cells = ncol(counts),
    has_bam = !is.na(bam) && file.exists(bam),
    has_barcode_file = file.exists(barcodes_file),
    stringsAsFactors = FALSE
  )

  rm(obj, counts, meta, meta_out)
  gc()
}

manifest <- bind_rows(manifest_rows) %>%
  arrange(.data$batch_type, .data$sample)

if (nrow(manifest) == 0) stop("No Numbat inputs were exported.")

bad <- manifest %>% filter(!.data$has_bam | !.data$has_barcode_file)
if (nrow(bad) > 0) {
  fwrite(bad, file.path(out_dir, "Auto_PDO_numbat_manifest_missing_inputs.csv"))
  stop(
    "Some samples are missing BAM or barcode inputs. See ",
    file.path(out_dir, "Auto_PDO_numbat_manifest_missing_inputs.csv")
  )
}

manifest_path <- file.path(out_dir, "Auto_PDO_numbat_manifest.csv")
fwrite(manifest, manifest_path)

summary_tbl <- manifest %>%
  count(.data$batch_type, .data$treatment, name = "n_samples") %>%
  arrange(.data$batch_type, .data$treatment)
fwrite(summary_tbl, file.path(out_dir, "Auto_PDO_numbat_manifest_summary.csv"))

message("Wrote Numbat manifest: ", file.path(out_root, manifest_path))
message("Samples exported: ", nrow(manifest))
