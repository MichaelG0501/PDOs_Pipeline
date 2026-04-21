####################
# Auto_3CA_pseudobulk_correlation_crossdata.R
#
# Compare pan-cancer 3CA metaprogram scores between scATLAS primary tumour
# and two PDO datasets using sample-level bulk-like inputs:
#   1. OAC PDO single-cell data pseudobulked by sample
#   2. OSCC/ESCC PDO bulk RNA-seq (GSE269447 tumor organoids only)
#
# IMPORTANT DOWNLOAD NOTE:
#   This script does NOT download OSCC bulk RNA-seq data.
#   Before running, download and extract GEO GSE269447 outside this repo to:
#   /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/
#   Expected files:
#     /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/GSE269447_RAW.tar
#     /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/raw_txt/GSM*_Tumor-Org_TPM.txt.gz
#   Example:
#     mkdir -p /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/raw_txt
#     curl -L -o /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/GSE269447_RAW.tar \
#       https://ftp.ncbi.nlm.nih.gov/geo/series/GSE269nnn/GSE269447/suppl/GSE269447_RAW.tar
#     tar -xf /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/GSE269447_RAW.tar \
#       -C /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/raw_txt
#
# Inputs:
#   PDOs_outs/PDOs_merged.rds
#   scRef_Pipeline/ref_outs/EAC_Ref_epi.rds
#   /rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv
#   /rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/raw_txt/GSM*_Tumor-Org_TPM.txt.gz
#
# Outputs:
#   PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata.pdf
#   PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata.png
#   PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata_summary.csv
#   PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata_mean_scores.csv
#   PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata_gene_overlap.csv
#   PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_scores_{pdo,scatlas,oscc}.csv
####################

library(Seurat)
library(UCell)
library(Matrix)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

out_dir <- "Auto_3CA_pseudobulk_correlation_crossdata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

oscc_dir <- "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/raw_txt"
three_ca_csv <- "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv"
excluded_pdo_samples <- "SUR843T3_PDO"
label_threshold <- 0.1

####################
# path helpers
####################
resolve_first_existing <- function(candidates, desc) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(desc, " not found in expected locations: ", paste(candidates, collapse = "; "))
  }
  existing[1]
}

load_3ca_gene_sets <- function(path) {
  mp_df <- read.csv(path, check.names = FALSE)
  mp_list <- as.list(mp_df)
  mp_list <- lapply(mp_list, function(x) unique(x[x != "" & !is.na(x)]))
  names(mp_list) <- make.names(sub("^MP", "3CA_mp_", names(mp_list)))
  mp_list
}

collapse_expression_to_symbols <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"

  feature_names <- rownames(expr_mat)
  ens_idx <- grepl("^ENSG", feature_names)
  collapsed_names <- feature_names

  if (any(ens_idx)) {
    ens_keys <- sub("\\..*$", "", feature_names[ens_idx])
    map_df <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = unique(ens_keys),
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    )
    map_df <- map_df[!is.na(map_df$SYMBOL) & map_df$SYMBOL != "", , drop = FALSE]
    map_df <- map_df[!duplicated(map_df$ENSEMBL), , drop = FALSE]
    symbol_map <- setNames(map_df$SYMBOL, map_df$ENSEMBL)
    mapped <- unname(symbol_map[ens_keys])
    use_original <- is.na(mapped) | mapped == ""
    collapsed_names[ens_idx] <- ifelse(use_original, feature_names[ens_idx], mapped)
  }

  keep_idx <- !is.na(collapsed_names) & collapsed_names != ""
  expr_mat <- expr_mat[keep_idx, , drop = FALSE]
  collapsed_names <- collapsed_names[keep_idx]

  rowsum(expr_mat, group = collapsed_names, reorder = FALSE)
}

make_sample_pseudobulk <- function(seurat_obj, exclude_samples = character(0)) {
  counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  sample_ids <- as.character(seurat_obj$orig.ident)
  keep_cells <- rep(TRUE, length(sample_ids))

  if (length(exclude_samples) > 0) {
    keep_cells <- !sample_ids %in% exclude_samples
    counts_mat <- counts_mat[, keep_cells, drop = FALSE]
    sample_ids <- sample_ids[keep_cells]
  }

  sample_factor <- factor(sample_ids, levels = unique(sample_ids))
  sample_mm <- sparse.model.matrix(~ 0 + sample_factor)
  colnames(sample_mm) <- levels(sample_factor)
  pseudobulk_counts <- counts_mat %*% sample_mm

  collapse_expression_to_symbols(pseudobulk_counts)
}

read_oscc_bulk_counts <- function(data_dir) {
  files <- sort(list.files(
    data_dir,
    pattern = "_Tumor-Org_TPM\\.txt\\.gz$",
    full.names = TRUE
  ))

  if (length(files) == 0) {
    stop("No tumor-organoid OSCC files found in: ", data_dir)
  }

  sample_vectors <- lapply(files, function(file_path) {
    dt <- fread(cmd = paste("gzip -dc", shQuote(file_path)), check.names = FALSE)
    count_col <- names(dt)[2]
    counts_vec <- dt[[count_col]]
    names(counts_vec) <- dt[[1]]
    counts_vec
  })

  all_genes <- unique(unlist(lapply(sample_vectors, names), use.names = FALSE))
  oscc_counts <- matrix(
    0,
    nrow = length(all_genes),
    ncol = length(sample_vectors),
    dimnames = list(all_genes, NULL)
  )

  for (idx in seq_along(sample_vectors)) {
    oscc_counts[names(sample_vectors[[idx]]), idx] <- sample_vectors[[idx]]
  }

  colnames(oscc_counts) <- vapply(files, function(file_path) {
    sample_id <- sub("_TPM\\.txt\\.gz$", "", basename(file_path))
    sample_id <- sub("^GSM\\d+_", "ESCCO", sample_id)
    sub("_Tumor-Org$", "_Tumor_Org", sample_id)
  }, character(1))

  collapse_expression_to_symbols(oscc_counts)
}

prepare_gene_sets <- function(expr_mat, gene_sets, dataset_name, min_genes = 5) {
  overlap_tbl <- lapply(names(gene_sets), function(mp_name) {
    genes_present <- intersect(gene_sets[[mp_name]], rownames(expr_mat))
    data.frame(
      dataset = dataset_name,
      mp = mp_name,
      genes_total = length(unique(gene_sets[[mp_name]])),
      genes_present = length(unique(genes_present)),
      prop_present = length(unique(genes_present)) / length(unique(gene_sets[[mp_name]])),
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  keep_mps <- overlap_tbl$mp[overlap_tbl$genes_present >= min_genes]
  gene_sets_trimmed <- lapply(gene_sets[keep_mps], function(gs) intersect(unique(gs), rownames(expr_mat)))

  list(
    gene_sets = gene_sets_trimmed,
    overlap = overlap_tbl
  )
}

clean_3ca_label <- function(x) {
  x <- sub("^X3CA_mp_", "", x)
  gsub("\\.", " ", x)
}

score_3ca_sets <- function(expr_mat, gene_sets) {
  as.data.frame(
    ScoreSignatures_UCell(
      matrix = expr_mat,
      features = gene_sets,
      name = "",
      ncores = 1
    )
  )
}

build_comparison_df <- function(target_mean, ref_mean, comparison_name) {
  common_mps <- intersect(names(target_mean), names(ref_mean))

  data.frame(
    MP = common_mps,
    scATLAS_score = ref_mean[common_mps],
    target_score = target_mean[common_mps],
    comparison = comparison_name,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Label = ifelse(
        scATLAS_score >= label_threshold | target_score >= label_threshold,
        clean_3ca_label(MP),
        NA
      ),
      Status = ifelse(
        scATLAS_score < label_threshold & target_score < label_threshold,
        "Low",
        "Highlighted"
      )
    )
}

format_p_value <- function(p_value) {
  out <- rep("NA", length(p_value))
  non_na_idx <- !is.na(p_value)
  sci_idx <- non_na_idx & p_value < 1e-4
  dec_idx <- non_na_idx & !sci_idx
  out[sci_idx] <- format(p_value[sci_idx], scientific = TRUE, digits = 2)
  out[dec_idx] <- sprintf("%.4f", p_value[dec_idx])
  out
}

####################
# Load data
####################
if (!file.exists(three_ca_csv)) {
  stop("3CA MP gene list file not found: ", three_ca_csv)
}

if (!dir.exists(oscc_dir)) {
  stop("OSCC data directory not found: ", oscc_dir)
}

message("Loading 3CA metaprogram gene sets...")
three_ca_gene_sets <- load_3ca_gene_sets(three_ca_csv)

####################
# resolve input paths with fallbacks
####################
pdo_rds_path <- resolve_first_existing(
  candidates = c(
    "PDOs_merged.rds",
    "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/PDOs_merged.rds",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/PDOs_merged.rds",
    "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/PDOs_merged.rds"
  ),
  desc = "PDOs_merged.rds"
)

scatlas_rds <- resolve_first_existing(
  candidates = c(
    "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds",
    "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds"
  ),
  desc = "EAC_Ref_epi.rds"
)

message("Loading OAC PDO Seurat object from: ", pdo_rds_path)
pdo_obj <- readRDS(pdo_rds_path)

message("Loading scATLAS Seurat object from: ", scatlas_rds)
scatlas_obj <- readRDS(scatlas_rds)

####################
# Build sample-level bulk-like matrices
####################
message("Building OAC PDO pseudobulk counts...")
pdo_pseudobulk <- make_sample_pseudobulk(
  pdo_obj,
  exclude_samples = excluded_pdo_samples
)

message("Building scATLAS pseudobulk counts...")
scatlas_pseudobulk <- make_sample_pseudobulk(scatlas_obj)

message("Loading OSCC bulk tumor organoid counts...")
oscc_bulk_counts <- read_oscc_bulk_counts(oscc_dir)

####################
# Restrict to comparable gene sets
####################
pdo_gene_sets <- prepare_gene_sets(pdo_pseudobulk, three_ca_gene_sets, "OAC PDO")
scatlas_gene_sets <- prepare_gene_sets(scatlas_pseudobulk, three_ca_gene_sets, "scATLAS")
oscc_gene_sets <- prepare_gene_sets(oscc_bulk_counts, three_ca_gene_sets, "OSCC PDO")

common_mps <- Reduce(intersect, list(
  names(pdo_gene_sets$gene_sets),
  names(scatlas_gene_sets$gene_sets),
  names(oscc_gene_sets$gene_sets)
))

if (length(common_mps) == 0) {
  stop("No shared 3CA MPs passed the minimum overlap filter across all datasets.")
}

pdo_gene_sets$gene_sets <- pdo_gene_sets$gene_sets[common_mps]
scatlas_gene_sets$gene_sets <- scatlas_gene_sets$gene_sets[common_mps]
oscc_gene_sets$gene_sets <- oscc_gene_sets$gene_sets[common_mps]

gene_overlap_df <- bind_rows(
  pdo_gene_sets$overlap,
  scatlas_gene_sets$overlap,
  oscc_gene_sets$overlap
) %>%
  filter(mp %in% common_mps) %>%
  arrange(factor(dataset, levels = c("scATLAS", "OAC PDO", "OSCC PDO")), mp)

write.csv(
  gene_overlap_df,
  file.path(out_dir, "Auto_3CA_pseudobulk_correlation_crossdata_gene_overlap.csv"),
  row.names = FALSE
)

####################
# Score shared 3CA MPs in each dataset
####################
message("Scoring OAC PDO pseudobulk samples...")
pdo_scores <- score_3ca_sets(pdo_pseudobulk, pdo_gene_sets$gene_sets)

message("Scoring scATLAS pseudobulk samples...")
scatlas_scores <- score_3ca_sets(scatlas_pseudobulk, scatlas_gene_sets$gene_sets)

message("Scoring OSCC PDO bulk samples...")
oscc_scores <- score_3ca_sets(oscc_bulk_counts, oscc_gene_sets$gene_sets)

write.csv(
  tibble::rownames_to_column(pdo_scores, "sample_id"),
  file.path(out_dir, "Auto_3CA_pseudobulk_scores_pdo.csv"),
  row.names = FALSE
)
write.csv(
  tibble::rownames_to_column(scatlas_scores, "sample_id"),
  file.path(out_dir, "Auto_3CA_pseudobulk_scores_scatlas.csv"),
  row.names = FALSE
)
write.csv(
  tibble::rownames_to_column(oscc_scores, "sample_id"),
  file.path(out_dir, "Auto_3CA_pseudobulk_scores_oscc.csv"),
  row.names = FALSE
)

####################
# Dataset-level mean scores and correlations
####################
pdo_mean_scores <- colMeans(pdo_scores[, common_mps, drop = FALSE], na.rm = TRUE)
scatlas_mean_scores <- colMeans(scatlas_scores[, common_mps, drop = FALSE], na.rm = TRUE)
oscc_mean_scores <- colMeans(oscc_scores[, common_mps, drop = FALSE], na.rm = TRUE)

mean_score_df <- data.frame(
  MP = common_mps,
  scATLAS = scatlas_mean_scores[common_mps],
  OAC_PDO = pdo_mean_scores[common_mps],
  OSCC_PDO = oscc_mean_scores[common_mps],
  stringsAsFactors = FALSE
) %>%
  mutate(MP_label = clean_3ca_label(MP))

write.csv(
  mean_score_df,
  file.path(out_dir, "Auto_3CA_pseudobulk_correlation_crossdata_mean_scores.csv"),
  row.names = FALSE
)

pdo_vs_sc <- build_comparison_df(
  target_mean = pdo_mean_scores,
  ref_mean = scatlas_mean_scores,
  comparison_name = "OAC PDO vs scATLAS"
)

oscc_vs_sc <- build_comparison_df(
  target_mean = oscc_mean_scores,
  ref_mean = scatlas_mean_scores,
  comparison_name = "OSCC PDO vs scATLAS"
)

####################
# panel labels with target sample counts
####################
comparison_levels <- c("OAC PDO vs scATLAS", "OSCC PDO vs scATLAS")
comparison_panels <- c(
  "OAC PDO vs scATLAS" = paste0("OAC PDO (n=", nrow(pdo_scores), ") vs scATLAS"),
  "OSCC PDO vs scATLAS" = paste0("OSCC PDO (n=", nrow(oscc_scores), ") vs scATLAS")
)

plot_df <- bind_rows(pdo_vs_sc, oscc_vs_sc) %>%
  mutate(
    comparison = factor(comparison, levels = comparison_levels),
    comparison_panel = factor(
      unname(comparison_panels[as.character(comparison)]),
      levels = unname(comparison_panels[comparison_levels])
    )
  )

summary_df <- plot_df %>%
  group_by(comparison) %>%
  group_modify(~ {
    cor_out <- cor.test(.x$scATLAS_score, .x$target_score, method = "spearman", exact = FALSE)
    tibble(
      rho = unname(cor_out$estimate),
      p_value = cor_out$p.value,
      n_mps = nrow(.x),
      highlight_mps = sum(.x$Status == "Highlighted"),
      max_scATLAS = max(.x$scATLAS_score, na.rm = TRUE),
      max_target = max(.x$target_score, na.rm = TRUE)
    )
  }) %>%
  ungroup() %>%
  mutate(
    target_dataset = ifelse(comparison == "OAC PDO vs scATLAS", "OAC PDO", "OSCC PDO"),
    target_sample_n = ifelse(comparison == "OAC PDO vs scATLAS", nrow(pdo_scores), nrow(oscc_scores)),
    scATLAS_sample_n = nrow(scatlas_scores)
  )

write.csv(
  summary_df,
  file.path(out_dir, "Auto_3CA_pseudobulk_correlation_crossdata_summary.csv"),
  row.names = FALSE
)

max_limit <- max(c(plot_df$scATLAS_score, plot_df$target_score), na.rm = TRUE) * 1.05
annotation_df <- summary_df %>%
  mutate(
    comparison_panel = factor(
      unname(comparison_panels[as.character(comparison)]),
      levels = unname(comparison_panels[comparison_levels])
    ),
    x = max_limit * 0.03,
    y = max_limit * 0.98,
    label = paste0(
      "Spearman rho = ", sprintf("%.2f", rho),
      "\nP = ", format_p_value(p_value)
    )
  )

####################
# Plot side-by-side correlations
####################
scatter_plot <- ggplot(plot_df, aes(x = scATLAS_score, y = target_score)) +
  geom_vline(
    xintercept = label_threshold,
    linetype = "dotted",
    color = "black",
    linewidth = 0.4,
    alpha = 0.5
  ) +
  geom_hline(
    yintercept = label_threshold,
    linetype = "dotted",
    color = "black",
    linewidth = 0.4,
    alpha = 0.5
  ) +
  geom_point(aes(color = Status), size = 2.8, alpha = 0.8) +
  geom_text_repel(
    aes(label = Label),
    size = 2.8,
    max.overlaps = 40,
    min.segment.length = 0,
    na.rm = TRUE
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#B2182B",
    linetype = "dashed",
    fill = "#B2182B",
    alpha = 0.12
  ) +
  geom_text(
    data = annotation_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.4
  ) +
  scale_color_manual(values = c("Low" = "grey65", "Highlighted" = "black")) +
  facet_wrap(~comparison_panel, nrow = 1) +
  coord_cartesian(xlim = c(0, max_limit), ylim = c(0, max_limit), expand = FALSE) +
  labs(
    title = "Pan-cancer 3CA MPs scATLAS vs PDOs correlation",
    subtitle = "OAC PDO and scATLAS pseudobulk, OSCC PDO from GSE269447 bulk data",
    x = "scATALS mean UCell score",
    y = "PDO mean UCell score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(face = "bold"),
    panel.spacing = grid::unit(1.2, "cm"),
    plot.margin = margin(10, 24, 10, 24)
  )

ggsave(
  file.path(out_dir, "Auto_3CA_pseudobulk_correlation_crossdata.pdf"),
  scatter_plot,
  width = 18,
  height = 7.5,
  useDingbats = FALSE
)

ggsave(
  file.path(out_dir, "Auto_3CA_pseudobulk_correlation_crossdata.png"),
  scatter_plot,
  width = 18,
  height = 7.5,
  dpi = 300
)

message("Completed Auto_3CA_pseudobulk_correlation_crossdata.R")
