####################
# Analysis registry:
#   Status: active terminal annotation export
#   Script: analysis/metaprograms/centred/Auto_06_centred_refined_mp_annotation_excel_export.R
#   Methodology: none; direct table export with no inferential or fitted thresholds
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Exports the canonical merged centred-refined MP gene lists in finalized
#     biological display order, with the current standardized MP descriptions.
#   Inputs:
#     - PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     - PDOs_outs/centred_mp_refinement/optimal_nMP.rds
#     - PDOs_outs/centred_mp_refinement/geneNMF_metaprograms_nMP_<optimal>.rds
#   Outputs:
#     - PDOs_outs/centred_mp_refinement/tables/merged_refined_MP_genes_summary.xlsx
#   Downstream use: none; terminal presentation/review workbook.
#   Cache/replot behavior: direct lightweight rebuild; no cache.
#   Run command:
#     Rscript analysis/metaprograms/centred/Auto_06_centred_refined_mp_annotation_excel_export.R
#   Conda env: dmtcp
####################
library(openxlsx)
library(dplyr)

live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
setwd(live_base)
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")

# Load inputs
merged_mp_genes <- readRDS("centred_mp_refinement/merged_refined_mp_genes.rds")
optimal_nMP <- readRDS("centred_mp_refinement/optimal_nMP.rds")
geneNMF.metaprograms <- readRDS(paste0("centred_mp_refinement/geneNMF_metaprograms_nMP_", optimal_nMP, ".rds"))

mp_desc_map <- setNames(
  paste(names(PDO_MP_DESCRIPTIONS), PDO_MP_DESCRIPTIONS, sep = "_"),
  names(PDO_MP_DESCRIPTIONS)
)

# Strip prefix for description
clean_desc <- function(mp_name) {
  if (mp_name %in% names(mp_desc_map)) {
    desc <- mp_desc_map[[mp_name]]
    sub("^MP[0-9]+[a-z\\+]*_", "", desc)
  } else {
    paste0(mp_name, "_unknown")
  }
}

get_desc <- function(mp_name) {
  clean_desc(mp_name)
}

merged_mps_ordered <- names(merged_mp_genes)

# Hard-coded user-specified ordering for exact MP order
user_mp_order <- c(PDO_CELL_CYCLE_MPS, unlist(PDO_MP_STATE_GROUPS, use.names = FALSE))
page1_order <- user_mp_order[user_mp_order %in% merged_mps_ordered]

parent_id <- function(x) {
  sub("\\+$", "", sub("[a-z]$", "", x))
}

# Logic for page 2 (Grouped_By_Parent)
# For by parent page, there should be small gap between distinct mps
# We can use the original tree order if possible, or just the sorted parent IDs
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order)]

all_parents <- parent_id(merged_mps_ordered)
unique_parents <- unique(all_parents)
unique_parent_ints <- as.integer(gsub("\\D", "", unique_parents))
ordered_unique_parents <- unique_parents[order(match(unique_parent_ints, mp_tree_order))]

page2_order <- character(0)
for (i in seq_along(ordered_unique_parents)) {
  p <- ordered_unique_parents[i]
  feats <- merged_mps_ordered[parent_id(merged_mps_ordered) == p]
  is_main <- !grepl("[a-z]$", feats)
  main_feat <- feats[is_main]
  if (length(main_feat) == 0) {
    main_feat <- p
  }
  sub_feats <- sort(feats[!is_main])
  page2_order <- c(page2_order, main_feat, sub_feats)
  if (i < length(ordered_unique_parents)) {
    page2_order <- c(page2_order, "GAP")
  }
}

build_mp_matrix <- function(mp_names_vec) {
  if (length(mp_names_vec) == 0) return(NULL)
  
  get_genes <- function(mp) {
    if (mp == "GAP") return(character(0))
    res <- merged_mp_genes[[mp]]
    if (is.null(res)) {
      nm <- gsub("\\+", "", mp)
      if (!is.null(geneNMF.metaprograms$metaprograms.genes[[nm]])) {
        res <- geneNMF.metaprograms$metaprograms.genes[[nm]]
      }
    }
    return(res)
  }
  
  max_g <- max(sapply(mp_names_vec, function(x) length(get_genes(x))))
  n_mp <- length(mp_names_vec)
  
  n_rows <- max_g + 2
  
  mat <- matrix(NA_character_, nrow = n_rows, ncol = n_mp)
  for (i in seq_along(mp_names_vec)) {
    mp <- mp_names_vec[i]
    if (mp == "GAP") {
      mat[1, i] <- ""
      mat[2, i] <- ""
    } else {
      mat[1, i] <- mp
      mat[2, i] <- get_desc(mp)
      genes <- get_genes(mp)
      if (length(genes) > 0) {
        mat[3:(length(genes)+2), i] <- genes
      }
    }
  }
  
  return(as.data.frame(mat, stringsAsFactors = FALSE))
}

# Create DataFrames
df_p1 <- build_mp_matrix(page1_order)
df_p2 <- build_mp_matrix(page2_order)

wb <- createWorkbook()

mp_name_style <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
desc_style <- createStyle(fgFill = "#F2F2F2")

add_sheet <- function(wb, sheet_name, df, order_vec) {
  addWorksheet(wb, sheet_name)
  sheet_idx <- length(names(wb))
  
  if (!is.null(df)) {
    # Write data starting from row 1 (no RETAINED/REMOVED header)
    writeData(wb, sheet = sheet_idx, x = df, startCol = 1, startRow = 1, colNames = FALSE)
    
    # Set styles and widths per column
    for (i in seq_along(order_vec)) {
      if (order_vec[i] == "GAP") {
        setColWidths(wb, sheet_idx, cols = i, widths = 3)
      } else {
        setColWidths(wb, sheet_idx, cols = i, widths = 25)
        addStyle(wb, sheet = sheet_idx, mp_name_style, rows = 1, cols = i, gridExpand = TRUE)
        addStyle(wb, sheet = sheet_idx, desc_style, rows = 2, cols = i, gridExpand = TRUE)
      }
    }
  }
}

add_sheet(wb, "Split_Separated", df_p1, page1_order)
add_sheet(wb, "Grouped_By_Parent", df_p2, page2_order)

out_dir <- "centred_mp_refinement/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(out_dir, "merged_refined_MP_genes_summary.xlsx")
saveWorkbook(wb, output_path, overwrite = TRUE)

message("Saved: ", output_path)
