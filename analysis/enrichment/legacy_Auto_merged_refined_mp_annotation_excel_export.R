####################
# Analysis registry (authoritative override):
#   Status: legacy; retained for provenance, no current downstream use
#   Script: analysis/enrichment/legacy_Auto_merged_refined_mp_annotation_excel_export.R
#   Methodology: historical method only; see analysis/ANALYSIS_MAP.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Preserves a superseded implementation or analysis tied to superseded
#     inputs. Do not use its outputs as current centred-MP/state inputs. The
#     original historical inputs, outputs, and method notes remain below.
####################

####################
# Analysis registry:
#   Status: active
#   Script: analysis/enrichment/Auto_merged_refined_mp_annotation_excel_export.R
#   Map: analysis/ANALYSIS_MAP.md
#   Description: 
#     Extract merged refined MP genes and create Excel summary.
#     MP description is currently just left as the raw name.
#     MP order strictly follows the unsupervised clustering order.
#   Inputs:
#     - live: PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#     - ephemeral: centred_mp_refinement/intermediate/merged_refined_mp_correlation_matrices.rds
#   Outputs (live: PDOs_outs/centred_mp_refinement/tables/):
#     merged_refined_MP_genes_summary_ordered.xlsx
####################

library(openxlsx)
library(dplyr)

project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
setwd(file.path(project_dir, "PDOs_outs"))

outdir_live <- "centred_mp_refinement"

# Load inputs
merged_mp_genes <- readRDS(file.path(outdir_live, "merged_refined_mp_genes.rds"))
cached_cor <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/centred_mp_refinement/intermediate/merged_refined_mp_correlation_matrices.rds"
cor_matrices <- readRDS(cached_cor)
mean_rho <- cor_matrices$mean_rho

library(ComplexHeatmap)

# Extract exact unsupervised clustering order (same as the plotting script)
# Need to use Heatmap to apply dendrogram reordering just like the plot does
ht_cor_unsup <- Heatmap(mean_rho, cluster_rows = TRUE, cluster_columns = TRUE)
pdf(NULL)
hm_drawn <- draw(ht_cor_unsup)
dev.off()
final_col_order <- colnames(mean_rho)[column_order(hm_drawn)]

# Ensure we only use available MPs
page1_order <- final_col_order[final_col_order %in% names(merged_mp_genes)]

get_desc <- function(mp_name) {
  # Currently leaving as raw name
  return(mp_name)
}

build_mp_matrix <- function(mp_names_vec) {
  if (length(mp_names_vec) == 0) return(NULL)
  
  get_genes <- function(mp) {
    if (mp == "GAP") return(character(0))
    res <- merged_mp_genes[[mp]]
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

# Create DataFrame using unsupervised clustering order
df_p1 <- build_mp_matrix(page1_order)

wb <- createWorkbook()

mp_name_style <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
desc_style <- createStyle(fgFill = "#F2F2F2")

add_sheet <- function(wb, sheet_name, df, order_vec) {
  addWorksheet(wb, sheet_name)
  sheet_idx <- length(names(wb))
  
  if (!is.null(df)) {
    writeData(wb, sheet = sheet_idx, x = df, startCol = 1, startRow = 1, colNames = FALSE)
    
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

add_sheet(wb, "Unsupervised_Clustered", df_p1, page1_order)

output_path <- file.path(outdir_live, "tables", "merged_refined_MP_genes_summary_ordered.xlsx")
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
saveWorkbook(wb, output_path, overwrite = TRUE)

message("Saved: ", output_path)
