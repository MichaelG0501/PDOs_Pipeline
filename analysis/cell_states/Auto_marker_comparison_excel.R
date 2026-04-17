####################
# Auto_marker_comparison_excel.R
#
# Cross-dataset comparison of state marker genes (scATLAS vs PDO).
# Reads saved ranked-marker and specificity data from both pipelines;
# produces one Excel sheet per shared state showing:
#   - Gene name
#   - scATLAS metrics (sample recurrence, study recurrence, effect size,
#     specificity gap, ranking score, ranking percent score)
#   - PDO metrics (sample recurrence, batch recurrence, effect size,
#     specificity gap, ranking score, ranking percent score)
#   - Combined ranking percent score
#   - Expression in the target state and other states for both datasets
#
# Inputs:
#   scRef: ref_outs/Auto_six_state_markers/Auto_six_state_markers_ranked.csv
#   scRef: ref_outs/Auto_six_state_markers/cache/state_specificity.rds
#   PDO:   PDOs_outs/Auto_five_state_markers/Auto_five_state_markers_ranked.csv
#   PDO:   PDOs_outs/Auto_five_state_markers/cache/state_specificity.rds
#
# Output:
#   PDOs_outs/Auto_five_state_markers/Auto_scATLAS_PDO_marker_comparison.xlsx
####################

####################
# libraries
####################
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(data.table)
})

####################
# paths
####################
scref_dir <- "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_six_state_markers"
pdo_dir   <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_five_state_markers"

out_xlsx <- file.path(pdo_dir, "Auto_scATLAS_PDO_marker_comparison.xlsx")

####################
# shared state mapping (scATLAS 6 states -> PDO 5 states)
####################
shared_states <- list(
  "Classic Proliferative" = list(
    sc = "Classic Proliferative",
    pdo = "Classic Proliferative"
  ),
  "Basal Metaplasia" = list(
    sc = "Basal to Intestinal Metaplasia",
    pdo = "Basal to Intest. Meta"
  ),
  "Stress-adaptive" = list(
    sc = "Stress-adaptive",
    pdo = "Stress-adaptive"
  ),
  "SMG-like Metaplasia" = list(
    sc = "SMG-like Metaplasia",
    pdo = "SMG-like Metaplasia"
  ),
  "3CA EMT & Protein Mat." = list(
    sc = "3CA_EMT_and_Protein_maturation",
    pdo = "3CA_EMT_and_Protein_maturation"
  )
)

# State order within each dataset (for expression columns)
sc_state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "3CA_EMT_and_Protein_maturation"
)

pdo_state_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "3CA_EMT_and_Protein_maturation"
)

# State colors for Excel headers
state_hex <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "3CA EMT & Protein Mat." = "#377EB8"
)

####################
# load data
####################
message("Loading ranked markers...")
sc_ranked  <- fread(file.path(scref_dir, "Auto_six_state_markers_ranked.csv"))
pdo_ranked <- fread(file.path(pdo_dir,   "Auto_five_state_markers_ranked.csv"))

message("Loading specificity caches (for per-state expression)...")
sc_spec  <- readRDS(file.path(scref_dir, "cache", "state_specificity.rds"))
pdo_spec <- readRDS(file.path(pdo_dir,   "cache", "state_specificity.rds"))

####################
# build expression matrices from specificity data
####################
build_expr_mat <- function(spec_df, state_order_use) {
  spec_df %>%
    select(gene, state, state_median_expr) %>%
    filter(state %in% state_order_use) %>%
    tidyr::pivot_wider(names_from = state, values_from = state_median_expr) %>%
    as.data.frame()
}

sc_expr_wide  <- build_expr_mat(sc_spec, sc_state_order)
pdo_expr_wide <- build_expr_mat(pdo_spec, pdo_state_order)

# Short display names for expression columns
sc_expr_display <- c(
  "Classic Proliferative" = "ClassProlif",
  "Basal to Intestinal Metaplasia" = "BasalMeta",
  "Stress-adaptive" = "StressAdapt",
  "SMG-like Metaplasia" = "SMG-like",
  "Immune Infiltrating" = "ImmuneInfil",
  "3CA_EMT_and_Protein_maturation" = "3CA_EMT"
)

pdo_expr_display <- c(
  "Classic Proliferative" = "ClassProlif",
  "Basal to Intest. Meta" = "BasalMeta",
  "Stress-adaptive" = "StressAdapt",
  "SMG-like Metaplasia" = "SMG-like",
  "3CA_EMT_and_Protein_maturation" = "3CA_EMT"
)

####################
# helper: normalise ranking score to 0-100 percent
####################
pct_score <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(50, length(x)))
  100 * (x - rng[1]) / diff(rng)
}

####################
# workbook creation
####################
wb <- createWorkbook()

# Styles
header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  border = "Bottom",
  borderStyle = "medium",
  wrapText = TRUE
)

sc_header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fontColour = "#FFFFFF",
  fgFill = "#2C3E50",
  border = "Bottom",
  borderStyle = "medium",
  wrapText = TRUE
)

pdo_header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fontColour = "#FFFFFF",
  fgFill = "#8E44AD",
  border = "Bottom",
  borderStyle = "medium",
  wrapText = TRUE
)

combined_header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fontColour = "#FFFFFF",
  fgFill = "#27AE60",
  border = "Bottom",
  borderStyle = "medium",
  wrapText = TRUE
)

expr_sc_header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fontColour = "#FFFFFF",
  fgFill = "#34495E",
  border = "Bottom",
  borderStyle = "medium",
  wrapText = TRUE
)

expr_pdo_header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  fontColour = "#FFFFFF",
  fgFill = "#7D3C98",
  border = "Bottom",
  borderStyle = "medium",
  wrapText = TRUE
)

num3 <- createStyle(numFmt = "0.000")
num2 <- createStyle(numFmt = "0.00")
num1 <- createStyle(numFmt = "0.0")
pct_fmt <- createStyle(numFmt = "0.0%")
sep_style <- createStyle(fgFill = "#D5D8DC", border = "LeftRight", borderColour = "#95A5A6")
gene_style <- createStyle(textDecoration = "bold", fontName = "Consolas")

# Alternating row bands
even_row_style <- createStyle(fgFill = "#F7F9F9")

for (st_label in names(shared_states)) {
  sc_state_name  <- shared_states[[st_label]]$sc
  pdo_state_name <- shared_states[[st_label]]$pdo

  message("Processing state: ", st_label)

  # Filter ranked markers to this state + significant DGE hits
  sc_hits <- sc_ranked %>%
    filter(state == sc_state_name) %>%
    as.data.frame()

  pdo_hits <- pdo_ranked %>%
    filter(state == pdo_state_name) %>%
    as.data.frame()

  # These are already filtered to hit_sample_n > 0, best_state_match, specificity_gap > 0
  # So they ARE significant DGE hits by construction of the pipeline

  # Compute percent ranking score within each dataset for this state
  if (nrow(sc_hits) > 0) {
    sc_hits$ranking_pct <- pct_score(sc_hits$ranking_score)
  }
  if (nrow(pdo_hits) > 0) {
    pdo_hits$ranking_pct <- pct_score(pdo_hits$ranking_score)
  }

  # Union of genes across both datasets
  all_genes <- sort(unique(c(sc_hits$gene, pdo_hits$gene)))
  if (length(all_genes) == 0) {
    message("  -> No significant markers for state: ", st_label, ". Skipping.")
    next
  }

  # Build scATLAS columns
  sc_df <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
  if (nrow(sc_hits) > 0) {
    sc_merge <- sc_hits %>%
      select(
        gene,
        sc_SampleRecur = sample_recurrence,
        sc_StudyRecur = study_recurrence,
        sc_log2FC = median_log2FC_hit,
        sc_SpecGap = specificity_gap,
        sc_RankScore = ranking_score,
        sc_RankPct = ranking_pct
      )
    sc_df <- sc_df %>% left_join(sc_merge, by = c("Gene" = "gene"))
  } else {
    sc_df$sc_SampleRecur <- NA_real_
    sc_df$sc_StudyRecur <- NA_real_
    sc_df$sc_log2FC <- NA_real_
    sc_df$sc_SpecGap <- NA_real_
    sc_df$sc_RankScore <- NA_real_
    sc_df$sc_RankPct <- NA_real_
  }

  # Build PDO columns
  if (nrow(pdo_hits) > 0) {
    pdo_merge <- pdo_hits %>%
      select(
        gene,
        pdo_SampleRecur = sample_recurrence,
        pdo_BatchRecur = batch_recurrence,
        pdo_log2FC = median_log2FC_hit,
        pdo_SpecGap = specificity_gap,
        pdo_RankScore = ranking_score,
        pdo_RankPct = ranking_pct
      )
    sc_df <- sc_df %>% left_join(pdo_merge, by = c("Gene" = "gene"))
  } else {
    sc_df$pdo_SampleRecur <- NA_real_
    sc_df$pdo_BatchRecur <- NA_real_
    sc_df$pdo_log2FC <- NA_real_
    sc_df$pdo_SpecGap <- NA_real_
    sc_df$pdo_RankScore <- NA_real_
    sc_df$pdo_RankPct <- NA_real_
  }

  # Combined ranking: average of percent scores (NA -> take the non-NA one)
  sc_df$CombinedRankPct <- rowMeans(
    cbind(sc_df$sc_RankPct, sc_df$pdo_RankPct),
    na.rm = TRUE
  )

  # Empty separator column
  sc_df$Sep <- ""

  # Expression columns - scATLAS
  sc_expr_cols <- intersect(sc_state_order, colnames(sc_expr_wide))
  if (length(sc_expr_cols) > 0 && nrow(sc_expr_wide) > 0) {
    # Subset for target state first then others
    target_sc_col <- sc_state_name
    other_sc_cols <- setdiff(sc_expr_cols, target_sc_col)
    sc_expr_ordered <- c(target_sc_col, other_sc_cols)

    sc_expr_sub <- sc_expr_wide[match(all_genes, sc_expr_wide$gene), c("gene", sc_expr_ordered), drop = FALSE]
    sc_expr_sub$gene <- NULL
    colnames(sc_expr_sub) <- paste0("sc_", sc_expr_display[sc_expr_ordered])
  } else {
    sc_expr_sub <- data.frame(matrix(NA, nrow = length(all_genes), ncol = 0))
  }

  # Expression columns - PDO
  pdo_expr_cols <- intersect(pdo_state_order, colnames(pdo_expr_wide))
  if (length(pdo_expr_cols) > 0 && nrow(pdo_expr_wide) > 0) {
    target_pdo_col <- pdo_state_name
    other_pdo_cols <- setdiff(pdo_expr_cols, target_pdo_col)
    pdo_expr_ordered <- c(target_pdo_col, other_pdo_cols)

    pdo_expr_sub <- pdo_expr_wide[match(all_genes, pdo_expr_wide$gene), c("gene", pdo_expr_ordered), drop = FALSE]
    pdo_expr_sub$gene <- NULL
    colnames(pdo_expr_sub) <- paste0("pdo_", pdo_expr_display[pdo_expr_ordered])
  } else {
    pdo_expr_sub <- data.frame(matrix(NA, nrow = length(all_genes), ncol = 0))
  }

  # Final combined table
  out_df <- cbind(
    sc_df[, c("Gene",
              "sc_SampleRecur", "sc_StudyRecur", "sc_log2FC", "sc_SpecGap", "sc_RankScore",
              "sc_RankPct",
              "pdo_SampleRecur", "pdo_BatchRecur", "pdo_log2FC", "pdo_SpecGap", "pdo_RankScore",
              "pdo_RankPct",
              "CombinedRankPct",
              "Sep"), drop = FALSE],
    sc_expr_sub,
    pdo_expr_sub
  )

  # Sort by combined ranking score descending
  out_df <- out_df[order(-out_df$CombinedRankPct, na.last = TRUE), ]

  # Column display names for header row 2
  display_names <- c(
    "Gene",
    "Sample\nRecurrence", "Study\nRecurrence", "Median\nlog2FC", "Specificity\nGap", "Ranking\nScore",
    "Rank %",
    "Sample\nRecurrence", "Batch\nRecurrence", "Median\nlog2FC", "Specificity\nGap", "Ranking\nScore",
    "Rank %",
    "Combined\nRank %",
    ""
  )
  # Add expression column display names
  if (ncol(sc_expr_sub) > 0) {
    sc_expr_names <- colnames(sc_expr_sub)
    sc_expr_names <- gsub("^sc_", "", sc_expr_names)
    # Highlight target state
    sc_expr_names[1] <- paste0(sc_expr_names[1], "\n(TARGET)")
    display_names <- c(display_names, sc_expr_names)
  }
  if (ncol(pdo_expr_sub) > 0) {
    pdo_expr_names <- colnames(pdo_expr_sub)
    pdo_expr_names <- gsub("^pdo_", "", pdo_expr_names)
    pdo_expr_names[1] <- paste0(pdo_expr_names[1], "\n(TARGET)")
    display_names <- c(display_names, pdo_expr_names)
  }

  # Sheet name (Excel limit 31 chars)
  sheet_name <- substr(st_label, 1, 31)
  addWorksheet(wb, sheet_name)

  # Compute column ranges
  n_total_cols <- ncol(out_df)
  n_data_rows <- nrow(out_df)

  # Row 1: Super-headers (dataset labels)
  # Gene col = 1
  # scATLAS = cols 2-7
  # PDO = cols 8-13
  # Combined = col 14
  # Sep = col 15
  # scATLAS expr = cols 16 to (15 + ncol(sc_expr_sub))
  # PDO expr = next set

  sc_metric_start <- 2
  sc_metric_end <- 7
  pdo_metric_start <- 8
  pdo_metric_end <- 13
  combined_col <- 14
  sep_col <- 15
  sc_expr_start <- 16
  sc_expr_end <- 15 + ncol(sc_expr_sub)
  pdo_expr_start <- sc_expr_end + 1
  pdo_expr_end <- sc_expr_end + ncol(pdo_expr_sub)

  # Write super-header row 1 with specific placement
  writeData(wb, sheet_name, "Gene", startCol = 1, startRow = 1)
  
  # Place metric headers at the start of their groups (Row 1)
  if (sc_metric_end >= sc_metric_start)
    writeData(wb, sheet_name, "scATLAS", startCol = sc_metric_start, startRow = 1)
  if (pdo_metric_end >= pdo_metric_start)
    writeData(wb, sheet_name, "PDO", startCol = pdo_metric_start, startRow = 1)
    
  writeData(wb, sheet_name, "", startCol = combined_col, startRow = 1)
  writeData(wb, sheet_name, "", startCol = sep_col, startRow = 1)
  
  # Split Expression headers into two columns to minimize row height
  if (ncol(sc_expr_sub) > 0) {
    writeData(wb, sheet_name, "scATLAS", startCol = sc_expr_start, startRow = 1)
    writeData(wb, sheet_name, "Expression", startCol = sc_expr_start + 1, startRow = 1)
  }
  if (ncol(pdo_expr_sub) > 0) {
    writeData(wb, sheet_name, "PDO", startCol = pdo_expr_start, startRow = 1)
    writeData(wb, sheet_name, "Expression", startCol = pdo_expr_start + 1, startRow = 1)
  }

  # Merge Gene header (A1:A2) as requested
  mergeCells(wb, sheet_name, cols = 1, rows = 1:2)
  
  # (Removed mergeCells for dataset metric headers per user request)
  # mergeCells(wb, sheet_name, cols = sc_metric_start:sc_metric_end, rows = 1)
  # mergeCells(wb, sheet_name, cols = pdo_metric_start:pdo_metric_end, rows = 1)
  
  # (Removed mergeCells for expression headers per user request to avoid height issues)
  # if (ncol(sc_expr_sub) > 1)
  #   mergeCells(wb, sheet_name, cols = sc_expr_start:sc_expr_end, rows = 1)
  # if (ncol(pdo_expr_sub) > 1)
  #   mergeCells(wb, sheet_name, cols = pdo_expr_start:pdo_expr_end, rows = 1)

  # Style super-header row
  addStyle(wb, sheet_name, sc_header_style, rows = 1, cols = sc_metric_start:sc_metric_end, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, pdo_header_style, rows = 1, cols = pdo_metric_start:pdo_metric_end, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, combined_header_style, rows = 1, cols = combined_col, stack = TRUE)
  if (ncol(sc_expr_sub) > 0)
    addStyle(wb, sheet_name, expr_sc_header_style, rows = 1, cols = sc_expr_start:sc_expr_end, gridExpand = TRUE, stack = TRUE)
  if (ncol(pdo_expr_sub) > 0)
    addStyle(wb, sheet_name, expr_pdo_header_style, rows = 1, cols = pdo_expr_start:pdo_expr_end, gridExpand = TRUE, stack = TRUE)

  # Write sub-header row 2 (individual column names)
  for (j in seq_along(display_names)) {
    writeData(wb, sheet_name, display_names[j], startCol = j, startRow = 2)
  }
  addStyle(wb, sheet_name, header_style, rows = 2, cols = 1:n_total_cols, gridExpand = TRUE, stack = TRUE)
  # Color sub-headers to match super-headers
  addStyle(wb, sheet_name, sc_header_style, rows = 2, cols = sc_metric_start:sc_metric_end, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, pdo_header_style, rows = 2, cols = pdo_metric_start:pdo_metric_end, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, combined_header_style, rows = 2, cols = combined_col, stack = TRUE)
  if (ncol(sc_expr_sub) > 0)
    addStyle(wb, sheet_name, expr_sc_header_style, rows = 2, cols = sc_expr_start:sc_expr_end, gridExpand = TRUE, stack = TRUE)
  if (ncol(pdo_expr_sub) > 0)
    addStyle(wb, sheet_name, expr_pdo_header_style, rows = 2, cols = pdo_expr_start:pdo_expr_end, gridExpand = TRUE, stack = TRUE)

  # Write data starting row 3
  writeData(wb, sheet_name, out_df, startCol = 1, startRow = 3, colNames = FALSE)

  # Gene name styling
  addStyle(wb, sheet_name, gene_style, rows = 3:(n_data_rows + 2), cols = 1, stack = TRUE)

  # Number formatting
  # Recurrence columns (proportions -> percentage format)
  recur_cols <- c(sc_metric_start, sc_metric_start + 1, pdo_metric_start, pdo_metric_start + 1)
  addStyle(wb, sheet_name, pct_fmt, rows = 3:(n_data_rows + 2), cols = recur_cols, gridExpand = TRUE, stack = TRUE)

  # log2FC, specificity gap columns (3 decimals)
  fc_gap_cols <- c(sc_metric_start + 2, sc_metric_start + 3, pdo_metric_start + 2, pdo_metric_start + 3)
  addStyle(wb, sheet_name, num3, rows = 3:(n_data_rows + 2), cols = fc_gap_cols, gridExpand = TRUE, stack = TRUE)

  # Ranking score (2 decimals)
  rank_cols <- c(sc_metric_start + 4, pdo_metric_start + 4)
  addStyle(wb, sheet_name, num2, rows = 3:(n_data_rows + 2), cols = rank_cols, gridExpand = TRUE, stack = TRUE)

  # Rank percent (1 decimal)
  rankpct_cols <- c(sc_metric_end, pdo_metric_end, combined_col)
  addStyle(wb, sheet_name, num1, rows = 3:(n_data_rows + 2), cols = rankpct_cols, gridExpand = TRUE, stack = TRUE)

  # Expression columns (3 decimals)
  expr_cols <- c()
  if (ncol(sc_expr_sub) > 0) expr_cols <- c(expr_cols, sc_expr_start:sc_expr_end)
  if (ncol(pdo_expr_sub) > 0) expr_cols <- c(expr_cols, pdo_expr_start:pdo_expr_end)
  if (length(expr_cols) > 0) {
    addStyle(wb, sheet_name, num3, rows = 3:(n_data_rows + 2), cols = expr_cols, gridExpand = TRUE, stack = TRUE)
  }

  # Separator column styling
  addStyle(wb, sheet_name, sep_style, rows = 1:(n_data_rows + 2), cols = sep_col, stack = TRUE)

  # Alternating row bands for readability
  even_rows <- seq(4, n_data_rows + 2, by = 2)
  if (length(even_rows) > 0) {
    addStyle(wb, sheet_name, even_row_style, rows = even_rows, cols = 1:n_total_cols, gridExpand = TRUE, stack = TRUE)
  }

  # Colour scales for key metrics
  if (n_data_rows > 1) {
    data_rows <- 3:(n_data_rows + 2)

    # Recurrence: white -> dark green
    for (cc in recur_cols) {
      conditionalFormatting(
        wb, sheet_name, cols = cc, rows = data_rows,
        style = c("#FFFFFF", "#27AE60"),
        rule = c(0, 1),
        type = "colourScale"
      )
    }

    # log2FC: white -> dark orange
    for (cc in c(sc_metric_start + 2, pdo_metric_start + 2)) {
      conditionalFormatting(
        wb, sheet_name, cols = cc, rows = data_rows,
        style = c("#FFFFFF", "#E67E22"),
        rule = c(0, 2),
        type = "colourScale"
      )
    }

    # Specificity gap: white -> dark blue
    for (cc in c(sc_metric_start + 3, pdo_metric_start + 3)) {
      conditionalFormatting(
        wb, sheet_name, cols = cc, rows = data_rows,
        style = c("#FFFFFF", "#2980B9"),
        rule = c(0, 1),
        type = "colourScale"
      )
    }

    # Ranking score: white -> warm yellow
    for (cc in rank_cols) {
      conditionalFormatting(
        wb, sheet_name, cols = cc, rows = data_rows,
        style = c("#FFFFFF", "#F39C12"),
        rule = c(0, 3),
        type = "colourScale"
      )
    }

    # Rank % and Combined rank %: white -> deep teal
    for (cc in rankpct_cols) {
      conditionalFormatting(
        wb, sheet_name, cols = cc, rows = data_rows,
        style = c("#FFFFFF", "#1ABC9C"),
        rule = c(0, 100),
        type = "colourScale"
      )
    }

    # Expression columns: white -> professional red (#B22222)
    # Using 3-point rule to improve contrast in high-expression range (2.2 vs 2.8)
    if (length(expr_cols) > 0) {
      for (cc in expr_cols) {
        conditionalFormatting(
          wb, sheet_name, cols = cc, rows = data_rows,
          style = c("#FFFFFF", "#FB8A8A", "#B22222"),
          rule = c(0, 1.5, 3.5),
          type = "colourScale"
        )
      }
    }
  }

  # Column widths
  setColWidths(wb, sheet_name, cols = 1, widths = 14)
  setColWidths(wb, sheet_name, cols = 2:(sep_col - 1), widths = 11)
  setColWidths(wb, sheet_name, cols = sep_col, widths = 2)
  if (length(expr_cols) > 0) {
    setColWidths(wb, sheet_name, cols = expr_cols, widths = 11)
  }

  # Freeze panes: freeze first column and first two rows
  freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 2)

  # Left border at section boundaries
  section_borders <- c(sc_metric_start, pdo_metric_start, combined_col, sep_col)
  if (ncol(sc_expr_sub) > 0) section_borders <- c(section_borders, sc_expr_start)
  if (ncol(pdo_expr_sub) > 0) section_borders <- c(section_borders, pdo_expr_start)

  border_style <- createStyle(border = "Left", borderStyle = "medium", borderColour = "#2C3E50")
  for (bc in section_borders) {
    addStyle(wb, sheet_name, border_style, rows = 1:(n_data_rows + 2), cols = bc, stack = TRUE)
  }

  # Bold the state colour in the tab
  # (openxlsx doesn't support tab colours natively; skip)

  message("  -> ", st_label, ": ", n_data_rows, " genes written.")
}

####################
# save
####################
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved comparison Excel to: ", out_xlsx)
####################
####################
# Generate 3-page Heatmap PDF
####################
message("Generating 3-page heatmap PDF...")
library(ComplexHeatmap)
library(circlize)

row_zscore <- function(mat) {
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

row_zscore_na <- function(mat) {
  z <- t(apply(mat, 1, function(x) {
    if (all(is.na(x))) return(rep(0, length(x)))
    x_scaled <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
    x_scaled[!is.finite(x_scaled)] <- 0
    return(x_scaled)
  }))
  colnames(z) <- colnames(mat)
  z
}

# Heatmap 1: scAtlas 6 states
sc_hm_markers <- sc_ranked %>%
  filter(hit_sample_n > 0, best_state_match, specificity_gap > 0) %>%
  group_by(state) %>% slice_head(n = 5) %>% ungroup() %>%
  mutate(state = factor(state, levels = sc_state_order)) %>%
  arrange(state)

sc_hm_expr <- sc_spec %>%
  select(gene, state, state_median_expr) %>%
  tidyr::pivot_wider(names_from = state, values_from = state_median_expr) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()
sc_hm_expr <- sc_hm_expr[sc_hm_markers$gene, sc_state_order, drop = FALSE]

dup_counter_sc <- ave(seq_along(sc_hm_markers$gene), sc_hm_markers$gene, FUN = seq_along)
rnames_sc <- ifelse(dup_counter_sc == 1, sc_hm_markers$gene, paste0(sc_hm_markers$gene, "__", dup_counter_sc))
rownames(sc_hm_expr) <- rnames_sc
sc_hm_z <- row_zscore(sc_hm_expr)

sc_state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666"
)

row_ann_sc <- rowAnnotation(
  State = sc_hm_markers$state,
  Sample_support = anno_barplot(sc_hm_markers$hit_sample_pct * 100, gp = gpar(fill = "grey35", col = NA), border = FALSE, width = unit(15, "mm")),
  Study_support = anno_barplot(sc_hm_markers$hit_study_n, gp = gpar(fill = "grey55", col = NA), border = FALSE, width = unit(15, "mm")),
  col = list(State = sc_state_cols),
  show_annotation_name = c(State = FALSE, Sample_support = TRUE, Study_support = TRUE),
  annotation_label = c(State = "State", Sample_support = "Sample\nSupport", Study_support = "Study\nSupport"),
  annotation_name_side = "top", annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 9, lineheight = 0.9),
  simple_anno_size = unit(4, "mm")
)
top_ann_sc <- HeatmapAnnotation(
  Top_State = factor(sc_state_order, levels = sc_state_order),
  col = list(Top_State = sc_state_cols),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(4, "mm")
)

col_fun <- colorRamp2(c(-2.0, 0, 2.5), c("#1D4E89", "#F8F4EC", "#B22222"))

ht_sc <- Heatmap(
  sc_hm_z, name = "Row Z-score",
  top_annotation = top_ann_sc, left_annotation = row_ann_sc, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = sc_hm_markers$state, row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, border = TRUE
)


# Heatmap 2: PDO 5 states
pdo_hm_markers <- pdo_ranked %>%
  filter(hit_sample_n > 0, best_state_match, specificity_gap > 0) %>%
  group_by(state) %>% slice_head(n = 5) %>% ungroup() %>%
  mutate(state = factor(state, levels = pdo_state_order)) %>%
  arrange(state)

pdo_hm_expr <- pdo_spec %>%
  select(gene, state, state_median_expr) %>%
  tidyr::pivot_wider(names_from = state, values_from = state_median_expr) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()
pdo_hm_expr <- pdo_hm_expr[pdo_hm_markers$gene, pdo_state_order, drop = FALSE]

dup_counter_pdo <- ave(seq_along(pdo_hm_markers$gene), pdo_hm_markers$gene, FUN = seq_along)
rnames_pdo <- ifelse(dup_counter_pdo == 1, pdo_hm_markers$gene, paste0(pdo_hm_markers$gene, "__", dup_counter_pdo))
rownames(pdo_hm_expr) <- rnames_pdo
pdo_hm_z <- row_zscore(pdo_hm_expr)

pdo_state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

row_ann_pdo <- rowAnnotation(
  State = pdo_hm_markers$state,
  Sample_support = anno_barplot(pdo_hm_markers$hit_sample_pct * 100, gp = gpar(fill = "grey35", col = NA), border = FALSE, width = unit(15, "mm")),
  Batch_support = anno_barplot(pdo_hm_markers$hit_batch_n, gp = gpar(fill = "grey55", col = NA), border = FALSE, width = unit(15, "mm")),
  col = list(State = pdo_state_cols),
  show_annotation_name = c(State = FALSE, Sample_support = TRUE, Batch_support = TRUE),
  annotation_label = c(State = "State", Sample_support = "Sample\nSupport", Batch_support = "Batch\nSupport"),
  annotation_name_side = "top", annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 9, lineheight = 0.9),
  simple_anno_size = unit(4, "mm")
)
top_ann_pdo <- HeatmapAnnotation(
  Top_State = factor(pdo_state_order, levels = pdo_state_order),
  col = list(Top_State = pdo_state_cols),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(4, "mm")
)

ht_pdo <- Heatmap(
  pdo_hm_z, name = "Row Z-score",
  top_annotation = top_ann_pdo, left_annotation = row_ann_pdo, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = pdo_hm_markers$state, row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, border = TRUE
)


# Heatmap 3: Combined Heatmap
comb_hm_list <- list()
for (st_label in names(shared_states)) {
  sc_st <- shared_states[[st_label]]$sc
  pdo_st <- shared_states[[st_label]]$pdo
  
  sc_h <- sc_ranked %>% filter(state == sc_st, hit_sample_n > 0, specificity_gap > 0)
  pdo_h <- pdo_ranked %>% filter(state == pdo_st, hit_sample_n > 0, specificity_gap > 0)
  
  if(nrow(sc_h) > 0) sc_h$ranking_pct <- pct_score(sc_h$ranking_score)
  if(nrow(pdo_h) > 0) pdo_h$ranking_pct <- pct_score(pdo_h$ranking_score)
  
  shared_genes <- intersect(sc_h$gene, pdo_h$gene)
  if (length(shared_genes) > 0) {
    df <- data.frame(gene = shared_genes, state = sc_st, stringsAsFactors = FALSE)
    sc_match <- match(shared_genes, sc_h$gene)
    pdo_match <- match(shared_genes, pdo_h$gene)
    
    df$CombinedRankPct <- rowMeans(cbind(sc_h$ranking_pct[sc_match], pdo_h$ranking_pct[pdo_match]), na.rm = TRUE)
    df$sc_sample_pct <- sc_h$hit_sample_pct[sc_match]
    df$sc_study_n <- sc_h$hit_study_n[sc_match]
    df$pdo_sample_pct <- pdo_h$hit_sample_pct[pdo_match]
    df$pdo_batch_n <- pdo_h$hit_batch_n[pdo_match]
    
    df <- df %>% arrange(desc(CombinedRankPct)) %>% head(5)
    comb_hm_list[[sc_st]] <- df
  }
}

# Handle Immune Infiltrating separately
im_h <- sc_ranked %>% filter(state == "Immune Infiltrating", hit_sample_n > 0, specificity_gap > 0)
if (nrow(im_h) > 0) {
  im_h$ranking_pct <- pct_score(im_h$ranking_score)
  im_h <- im_h %>% arrange(desc(ranking_pct)) %>% head(5)
  df <- data.frame(
    gene = im_h$gene,
    state = "Immune Infiltrating",
    CombinedRankPct = im_h$ranking_pct,
    sc_sample_pct = im_h$hit_sample_pct,
    sc_study_n = im_h$hit_study_n,
    pdo_sample_pct = 0,
    pdo_batch_n = 0,
    stringsAsFactors = FALSE
  )
  comb_hm_list[["Immune Infiltrating"]] <- df
}

comb_markers <- bind_rows(comb_hm_list)
comb_markers$state <- factor(comb_markers$state, levels = sc_state_order)
comb_markers <- comb_markers %>% arrange(state)

sc_expr_all <- sc_spec %>%
  select(gene, state, state_median_expr) %>%
  tidyr::pivot_wider(names_from = state, values_from = state_median_expr) %>%
  tibble::column_to_rownames("gene") %>% as.matrix()

pdo_expr_all <- pdo_spec %>%
  select(gene, state, state_median_expr) %>%
  tidyr::pivot_wider(names_from = state, values_from = state_median_expr) %>%
  tibble::column_to_rownames("gene") %>% as.matrix()

comb_genes <- comb_markers$gene

mat_sc <- sc_expr_all[comb_genes, sc_state_order, drop=FALSE]
colnames(mat_sc) <- paste0("scAtlas::", sc_state_order)

mat_pdo <- matrix(NA, nrow=length(comb_genes), ncol=length(pdo_state_order))
rownames(mat_pdo) <- comb_genes
colnames(mat_pdo) <- paste0("PDO::", pdo_state_order)
genes_in_pdo <- comb_genes %in% rownames(pdo_expr_all)
mat_pdo[genes_in_pdo, ] <- pdo_expr_all[comb_genes[genes_in_pdo], pdo_state_order, drop=FALSE]

comb_expr <- cbind(mat_sc, mat_pdo)

mat_sc_z  <- row_zscore(mat_sc)
mat_pdo_z <- row_zscore_na(mat_pdo)
comb_z <- cbind(mat_sc_z, mat_pdo_z)

dup_counter_comb <- ave(seq_along(comb_genes), comb_genes, FUN = seq_along)
rnames_comb <- ifelse(dup_counter_comb == 1, comb_genes, paste0(comb_genes, "__", dup_counter_comb))
rownames(comb_z) <- rnames_comb

row_ann_comb <- rowAnnotation(
  State = comb_markers$state,
  sc_Sample = anno_barplot(comb_markers$sc_sample_pct * 100, gp = gpar(fill = "grey35", col = NA), border = FALSE, width = unit(15, "mm")),
  sc_Study = anno_barplot(comb_markers$sc_study_n, gp = gpar(fill = "grey55", col = NA), border = FALSE, width = unit(15, "mm")),
  pdo_Sample = anno_barplot(comb_markers$pdo_sample_pct * 100, gp = gpar(fill = "grey35", col = NA), border = FALSE, width = unit(15, "mm")),
  pdo_Batch = anno_barplot(comb_markers$pdo_batch_n, gp = gpar(fill = "grey55", col = NA), border = FALSE, width = unit(15, "mm")),
  col = list(State = sc_state_cols), 
  show_annotation_name = c(State=FALSE, sc_Sample=TRUE, sc_Study=TRUE, pdo_Sample=TRUE, pdo_Batch=TRUE),
  annotation_label = c(State="State", sc_Sample="scAtlas\nSample", sc_Study="scAtlas\nStudy", pdo_Sample="PDO\nSample", pdo_Batch="PDO\nBatch"),
  annotation_name_side = "top", annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize=8, lineheight=0.9),
  gap = unit(c(2, 2, 4, 2, 6), "mm"),
  simple_anno_size = unit(4, "mm")
)

col_split_factors <- factor(
  c(rep("scAtlas", length(sc_state_order)), rep("PDO", length(pdo_state_order))),
  levels = c("scAtlas", "PDO")
)

col_state_vec <- c(sc_state_order, pdo_state_order)
all_state_cols <- c(sc_state_cols, pdo_state_cols)

top_ann_comb <- HeatmapAnnotation(
  Dataset = col_split_factors,
  State = factor(col_state_vec, levels = unique(col_state_vec)),
  col = list(
    Dataset = c("scAtlas" = "#2C3E50", "PDO" = "#8E44AD"),
    State = all_state_cols
  ),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(4, "mm")
)

ht_comb <- Heatmap(
  comb_z, name = "Row Z-score",
  top_annotation = top_ann_comb, left_annotation = row_ann_comb, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = comb_markers$state, row_title_rot = 0,
  column_split = col_split_factors,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, border = TRUE
)


pdf(file.path(pdo_dir, "Auto_combined_marker_heatmaps.pdf"), width = 17, height = 12, useDingbats = FALSE)
draw(ht_sc)
draw(ht_pdo)
draw(ht_comb)
dev.off()
####################
# Output Top 5 Subset Excel (3 Sheets: scAtlas, PDO, Combined)
####################
# Build a separate top5 Excel with 3 sheets matching the SCENIC subset style
message("Generating Top 5 marker subset Excel (3 sheets)...")

# Ensure colors for all states
sc_state_cols_full <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666"
)

wb2 <- createWorkbook()

build_top5_markers_subset_sheet <- function(wb, sheet_name, mapping_df) {
  addWorksheet(wb, sheet_name)
  
  # Fix: Use a unified state order that mirrors scAtlas sequence for all datasets
  # (ensures 'Basal to Intest. Meta' is correctly ordered as the 2nd block)
  all_possible_states <- c(
    "Classic Proliferative", 
    "Basal to Intestinal Metaplasia", "Basal to Intest. Meta",
    "Stress-adaptive", 
    "SMG-like Metaplasia", 
    "Immune Infiltrating", 
    "3CA_EMT_and_Protein_maturation"
  )
  states_present <- intersect(all_possible_states, as.character(mapping_df$state))
  
  mapping_df$state <- factor(as.character(mapping_df$state), levels = states_present)
  mapping_df <- mapping_df %>% arrange(state)
  
  # Build final data frame with separators and grouped columns
  final_df <- data.frame()
  row_colors <- character()
  is_data_row <- logical()
  state_block_rows <- list() 
  
  current_excel_row <- 3
  
  for (i in seq_along(states_present)) {
    st <- states_present[i]
    st_genes <- mapping_df$gene[as.character(mapping_df$state) == st]
    if (length(st_genes) == 0) next
    
    # 1. Vertical State Label + Base Gene info
    st_data <- data.frame(
      StateLabel = rep(as.character(st), length(st_genes)),
      Gene = st_genes, 
      stringsAsFactors = FALSE
    )
    
    # Save block range
    block_start <- current_excel_row
    block_end <- current_excel_row + length(st_genes) - 1
    state_block_rows[[as.character(st)]] <- c(block_start, block_end)
    
    # 2. scAtlas Recurrence
    sc_r <- sc_ranked %>% filter(gene %in% st_genes) %>%
            select(gene, sc_SampleRecur = hit_sample_pct, sc_StudyRecur = hit_study_n)
    sc_r <- sc_r[!duplicated(sc_r$gene), ]
    st_data <- left_join(st_data, sc_r, by = c("Gene" = "gene"))
    
    # 3. PDO Recurrence
    pdo_r <- pdo_ranked %>% filter(gene %in% st_genes) %>%
             select(gene, pdo_SampleRecur = hit_sample_pct, pdo_BatchRecur = hit_batch_n)
    pdo_r <- pdo_r[!duplicated(pdo_r$gene), ]
    st_data <- left_join(st_data, pdo_r, by = c("Gene" = "gene"))
    
    # 4. Separator
    st_data$Sep <- ""
    
    # 5. scAtlas Expression Group
    sc_e <- sc_expr_wide[match(st_genes, sc_expr_wide$gene), sc_state_order, drop=FALSE]
    colnames(sc_e) <- paste0("sc_", sc_state_order)
    st_data <- cbind(st_data, sc_e)
    
    # 6. PDO Expression Group
    pdo_e <- pdo_expr_wide[match(st_genes, pdo_expr_wide$gene), pdo_state_order, drop=FALSE]
    colnames(pdo_e) <- paste0("pdo_", pdo_state_order)
    st_data <- cbind(st_data, pdo_e)
    
    final_df <- rbind(final_df, st_data)
    
    # Color
    st_col <- sc_state_cols_full[as.character(st)]
    if (is.na(st_col)) st_col <- "#000000"
    row_colors <- c(row_colors, rep(st_col, nrow(st_data)))
    is_data_row <- c(is_data_row, rep(TRUE, nrow(st_data)))
    
    current_excel_row <- block_end + 1
    
    # Empty Row
    if (i < length(states_present)) {
      empty_row <- as.data.frame(matrix(NA, nrow=1, ncol=ncol(st_data)))
      colnames(empty_row) <- colnames(st_data)
      empty_row$StateLabel <- ""
      empty_row$Gene <- ""
      empty_row$Sep <- ""
      final_df <- rbind(final_df, empty_row)
      row_colors <- c(row_colors, NA)
      is_data_row <- c(is_data_row, FALSE)
      current_excel_row <- current_excel_row + 1
    }
  }
  
  # Column indices
  gene_col <- 2
  sc_recur_cols <- 3:4
  pdo_recur_cols <- 5:6
  sep_col <- 7
  sc_expr_start <- 8
  sc_expr_end <- 7 + length(sc_state_order)
  pdo_expr_start <- sc_expr_end + 1
  pdo_expr_end <- sc_expr_end + length(pdo_state_order)
  
  # Row 1 Headers
  writeData(wb, sheet_name, "State", startCol=1, startRow=1)
  mergeCells(wb, sheet_name, cols=1, rows=1:2)
  writeData(wb, sheet_name, "Gene", startCol=2, startRow=1)
  mergeCells(wb, sheet_name, cols=2, rows=1:2)
  
  writeData(wb, sheet_name, "scAtlas Recurrence", startCol=sc_recur_cols[1], startRow=1)
  mergeCells(wb, sheet_name, cols=sc_recur_cols, rows=1)
  addStyle(wb, sheet_name, sc_header_style, rows=1, cols=sc_recur_cols, gridExpand=TRUE, stack=TRUE)
  
  writeData(wb, sheet_name, "PDO Recurrence", startCol=pdo_recur_cols[1], startRow=1)
  mergeCells(wb, sheet_name, cols=pdo_recur_cols, rows=1)
  addStyle(wb, sheet_name, pdo_header_style, rows=1, cols=pdo_recur_cols, gridExpand=TRUE, stack=TRUE)
  
  writeData(wb, sheet_name, "scAtlas Mean Expression", startCol=sc_expr_start, startRow=1)
  mergeCells(wb, sheet_name, cols=sc_expr_start:sc_expr_end, rows=1)
  addStyle(wb, sheet_name, expr_sc_header_style, rows=1, cols=sc_expr_start:sc_expr_end, gridExpand=TRUE, stack=TRUE)
  
  writeData(wb, sheet_name, "PDO Mean Expression", startCol=pdo_expr_start, startRow=1)
  mergeCells(wb, sheet_name, cols=pdo_expr_start:pdo_expr_end, rows=1)
  addStyle(wb, sheet_name, expr_pdo_header_style, rows=1, cols=pdo_expr_start:pdo_expr_end, gridExpand=TRUE, stack=TRUE)
  
  # Row 2 Headers
  headers2 <- c("State", "Gene", "Sample\nRecurrence", "Study\nRecurrence", "Sample\nRecurrence", "Batch\nRecurrence", "")
  headers2 <- c(headers2, sc_expr_display[sc_state_order], pdo_expr_display[pdo_state_order])
  writeData(wb, sheet_name, t(headers2), startCol=1, startRow=2, colNames=FALSE)
  addStyle(wb, sheet_name, header_style, rows=2, cols=1:ncol(final_df), gridExpand=TRUE, stack=TRUE)
  addStyle(wb, sheet_name, sc_header_style, rows = 2, cols = sc_recur_cols, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, pdo_header_style, rows = 2, cols = pdo_recur_cols, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, expr_sc_header_style, rows = 2, cols = sc_expr_start:sc_expr_end, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, expr_pdo_header_style, rows = 2, cols = pdo_expr_start:pdo_expr_end, gridExpand = TRUE, stack = TRUE)
  
  # Data
  writeData(wb, sheet_name, final_df, startCol = 1, startRow = 3, colNames = FALSE)
  
  # Styling for Column 1 (Merged Vertical State) and Column 2 (Gene color)
  for (st_name in names(state_block_rows)) {
    rows <- state_block_rows[[st_name]]
    mergeCells(wb, sheet_name, cols = 1, rows = rows[1]:rows[2])
    
    st_col_hex <- sc_state_cols_full[st_name]
    if (is.na(st_col_hex)) st_col_hex <- "#000000"
    
    vert_style <- createStyle(textRotation = 90, halign = "center", valign = "center", textDecoration = "bold", 
                              fontColour = "#FFFFFF", fgFill = st_col_hex, border = "TopBottomLeftRight")
    addStyle(wb, sheet_name, vert_style, rows = rows[1]:rows[2], cols = 1, stack = TRUE)
  }

  for (r in seq_len(nrow(final_df))) {
    if (is_data_row[r]) {
      st_col_hex <- row_colors[r]
      gene_st_style <- createStyle(textDecoration = "bold", fontName = "Consolas", fontColour = st_col_hex)
      addStyle(wb, sheet_name, gene_st_style, rows = r + 2, cols = 2, stack = TRUE)
    }
  }
  
  # Metric formatting
  addStyle(wb, sheet_name, pct_fmt, rows = 3:(nrow(final_df)+2), cols = c(3, 5), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, createStyle(numFmt = "0"), rows = 3:(nrow(final_df)+2), cols = c(4, 6), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, num3, rows = 3:(nrow(final_df)+2), cols = sc_expr_start:pdo_expr_end, gridExpand = TRUE, stack = TRUE)
  
  # Separator and borders
  addStyle(wb, sheet_name, sep_style, rows = 1:(nrow(final_df)+2), cols = sep_col, stack = TRUE)
  border_style_med <- createStyle(border = "Left", borderStyle = "medium", borderColour = "#2C3E50")
  addStyle(wb, sheet_name, border_style_med, rows = 1:(nrow(final_df)+2), cols = c(2, 3, 5, 7, sc_expr_start, pdo_expr_start), gridExpand = TRUE, stack = TRUE)
  
  # Conditional formatting targeting higher contrast (scAtlas)
  sc_expr_vals <- unlist(final_df[is_data_row, sc_expr_start:sc_expr_end])
  sc_expr_vals <- sc_expr_vals[is.finite(sc_expr_vals)]
  if (length(sc_expr_vals) > 0) {
    q_low <- quantile(sc_expr_vals, 0.05, na.rm=TRUE)
    q_mid <- quantile(sc_expr_vals, 0.40, na.rm=TRUE) # Pulled mid down for more red
    q_high <- quantile(sc_expr_vals, 0.85, na.rm=TRUE) # Lowered high to avoid outlier-lightening
    conditionalFormatting(wb, sheet_name, cols = sc_expr_start:sc_expr_end, rows = 3:(nrow(final_df)+2),
                          style = c("#FFFFFF", "#FB8A8A", "#B22222"), rule = c(q_low, q_mid, q_high), type = "colourScale")
  }
  
  # Conditional formatting (PDO)
  pdo_expr_vals <- unlist(final_df[is_data_row, pdo_expr_start:pdo_expr_end])
  pdo_expr_vals <- pdo_expr_vals[is.finite(pdo_expr_vals)]
  if (length(pdo_expr_vals) > 0) {
    q_low_p <- quantile(pdo_expr_vals, 0.05, na.rm=TRUE)
    q_mid_p <- quantile(pdo_expr_vals, 0.40, na.rm=TRUE)
    q_high_p <- quantile(pdo_expr_vals, 0.85, na.rm=TRUE)
    conditionalFormatting(wb, sheet_name, cols = pdo_expr_start:pdo_expr_end, rows = 3:(nrow(final_df)+2),
                          style = c("#FFFFFF", "#FB8A8A", "#B22222"), rule = c(q_low_p, q_mid_p, q_high_p), type = "colourScale")
  }
  
  # Layout
  setColWidths(wb, sheet_name, cols = 1, widths = 5)
  setColWidths(wb, sheet_name, cols = 2, widths = 20)
  setColWidths(wb, sheet_name, cols = 3:ncol(final_df), widths = 12)
  setColWidths(wb, sheet_name, cols = sep_col, widths = 3)
  freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 3)
}

# Generate 3 sheets
build_top5_markers_subset_sheet(wb2, "scAtlas Top 5", sc_hm_markers)
build_top5_markers_subset_sheet(wb2, "PDO Top 5", pdo_hm_markers)
build_top5_markers_subset_sheet(wb2, "Combined Top 5", comb_markers)

out_xlsx_top5 <- file.path(pdo_dir, "Auto_scATLAS_PDO_marker_top5.xlsx")
saveWorkbook(wb2, out_xlsx_top5, overwrite = TRUE)
message("Saved TOP 5 marker Excel (3 sheets) to: ", out_xlsx_top5)

