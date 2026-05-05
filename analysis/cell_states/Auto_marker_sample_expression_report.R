####################
# Auto_marker_sample_expression_report.R
#
# For the 3-per-state marker panels defined in Auto_marker_selection_simulation.R,
# produce:
#   1. An Excel workbook (styled like Auto_marker_comparison_excel.R) with four
#      data chunks — expression and pct-detected for SUR1090_Untreated_PDO and
#      SUR1072_Untreated_PDO — rows grouped by state.
#   2. A PDF dot plot similar to the simulation detection-by-state plot, but
#      showing two sample columns per state with raw expression (size = pct
#      detected, colour = value normalised across all 10 state×sample combos).
#
# Inputs:
#   PDOs_outs/PDOs_merged.rds
#   PDOs_outs/Auto_PDO_final_states.rds
#
# Outputs:
#   PDOs_outs/Auto_marker_selection_simulation/Auto_marker_sample_expression.xlsx
#   PDOs_outs/Auto_marker_selection_simulation/Auto_marker_sample_specificity_dotplot.pdf
####################

####################
# libraries
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(openxlsx)
  library(ggtext)
})

####################
# setup
####################
project_candidates <- c(
  Sys.getenv("PDO_PROJECT_DIR", unset = ""),
  "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline",
  "/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline",
  getwd()
)
project_candidates <- project_candidates[nzchar(project_candidates)]
project_dir <- project_candidates[dir.exists(project_candidates)][1]
if (is.na(project_dir)) stop("Could not locate the PDOs_Pipeline project directory.")

setwd(file.path(project_dir, "PDOs_outs"))
out_dir <- "Auto_marker_selection_simulation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

####################
# marker panels & constants
####################
state_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Total" = "black"
)

marker_panels <- list(
  "Classic Proliferative" = c("PCLAF", "STMN1", "FABP5"),
  "Basal to Intest. Meta" = c("MUC13", "HPGD", "TSPAN1"),
  "SMG-like Metaplasia"   = c("PROM1", "CHRM3", "PLCB4"),
  "Stress-adaptive"       = c("PPP1R15A", "ATF3", "SOX4"),
  "3CA_EMT_and_Protein_maturation" = c("CANX", "PDIA4", "NORAD")
)

samples_use <- c("SUR1090_Untreated_PDO", "SUR1072_Untreated_PDO")
sample_labels_short <- c(
  "SUR1090_Untreated_PDO" = "SUR1090 Untreated",
  "SUR1072_Untreated_PDO" = "SUR1072 Untreated"
)

####################
# load data
####################
message("Loading Seurat object...")
pdos <- readRDS("PDOs_merged.rds")
state_labels <- readRDS("Auto_PDO_final_states.rds")
DefaultAssay(pdos) <- "RNA"

common_cells <- intersect(colnames(pdos), names(state_labels))
state_labels <- state_labels[common_cells]

keep_cells <- common_cells[
  as.character(state_labels) %in% state_order &
    as.character(pdos$orig.ident[match(common_cells, colnames(pdos))]) %in% samples_use
]
cell_states <- factor(as.character(state_labels[keep_cells]), levels = state_order)
cell_samples <- as.character(pdos$orig.ident[match(keep_cells, colnames(pdos))])

all_marker_genes <- unique(unlist(marker_panels, use.names = FALSE))

get_assay_layer <- function(obj, assay = "RNA", layer = "counts") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
}

data_mat <- get_assay_layer(pdos, "RNA", "data")[all_marker_genes, keep_cells, drop = FALSE]
counts_mat <- get_assay_layer(pdos, "RNA", "counts")[all_marker_genes, keep_cells, drop = FALSE]

rm(pdos, state_labels); invisible(gc())

# compute per-state per-sample expression and pct detected
results <- list()
for (samp in samples_use) {
  # Per-state
  for (st in state_order) {
    idx <- which(cell_states == st & cell_samples == samp)
    n_cells <- length(idx)
    if (n_cells == 0) {
      expr_vals <- setNames(rep(NA_real_, length(all_marker_genes)), all_marker_genes)
      pct_vals  <- expr_vals
    } else {
      expr_vals <- Matrix::rowMeans(data_mat[, idx, drop = FALSE])
      pct_vals  <- Matrix::rowMeans(counts_mat[, idx, drop = FALSE] > 0)
    }
    results[[paste(st, samp, sep = "||")]] <- data.frame(
      gene = all_marker_genes,
      state = st,
      sample = samp,
      mean_expr = as.numeric(expr_vals[all_marker_genes]),
      pct_detected = as.numeric(pct_vals[all_marker_genes]),
      n_cells = n_cells,
      stringsAsFactors = FALSE
    )
  }
  # Combined (Total) - real overall expression (weighted average across all cells)
  idx_tot <- which(cell_samples == samp)
  n_cells_tot <- length(idx_tot)
  if (n_cells_tot == 0) {
    expr_vals_tot <- setNames(rep(NA_real_, length(all_marker_genes)), all_marker_genes)
    pct_vals_tot  <- expr_vals_tot
  } else {
    expr_vals_tot <- Matrix::rowMeans(data_mat[, idx_tot, drop = FALSE])
    pct_vals_tot  <- Matrix::rowMeans(counts_mat[, idx_tot, drop = FALSE] > 0)
  }
  results[[paste("Total", samp, sep = "||")]] <- data.frame(
    gene = all_marker_genes,
    state = "Total",
    sample = samp,
    mean_expr = as.numeric(expr_vals_tot[all_marker_genes]),
    pct_detected = as.numeric(pct_vals_tot[all_marker_genes]),
    n_cells = n_cells_tot,
    stringsAsFactors = FALSE
  )
}
res_df <- bind_rows(results)

####################
# ============ EXCEL WORKBOOK ============
####################
message("Building Excel workbook...")

# Build manifest ordered by state
manifest <- data.frame(
  panel = rep(names(marker_panels), lengths(marker_panels)),
  gene  = unlist(marker_panels, use.names = FALSE),
  stringsAsFactors = FALSE
)
manifest$panel <- factor(manifest$panel, levels = state_order)

wb <- createWorkbook()
addWorksheet(wb, "Marker Expression")
sn <- "Marker Expression"

# --- styles (matching Auto_marker_comparison_excel.R) ---
header_style <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
chunk1_hdr <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  fontColour = "#FFFFFF", fgFill = "#2C3E50",
  border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
chunk2_hdr <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  fontColour = "#FFFFFF", fgFill = "#8E44AD",
  border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
chunk3_hdr <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  fontColour = "#FFFFFF", fgFill = "#34495E",
  border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
chunk4_hdr <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  fontColour = "#FFFFFF", fgFill = "#7D3C98",
  border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
num3  <- createStyle(numFmt = "0.000")
pct_fmt <- createStyle(numFmt = "0.0%")
sep_style <- createStyle(fgFill = "#D5D8DC", border = "LeftRight", borderColour = "#95A5A6")
gene_style <- createStyle(textDecoration = "bold", fontName = "Consolas")
even_row_style <- createStyle(fgFill = "#F7F9F9")
border_left <- createStyle(border = "Left", borderStyle = "medium", borderColour = "#2C3E50")

# Pivot expression/pct for each sample × each state
state_order_ext <- c(state_order, "Total")
state_short <- c(
  "Classic Proliferative" = "ClassProlif",
  "Basal to Intest. Meta" = "BasalMeta",
  "SMG-like Metaplasia" = "SMG-like",
  "Stress-adaptive" = "StressAdapt",
  "3CA_EMT_and_Protein_maturation" = "3CA_EMT",
  "Total" = "Total"
)

# For each gene, get expr/pct per state for each sample
build_chunk <- function(samp, metric) {
  sub <- res_df %>% filter(sample == samp) %>%
    select(gene, state, value = !!sym(metric)) %>%
    mutate(state = factor(state, levels = state_order_ext)) %>%
    pivot_wider(names_from = state, values_from = value) %>%
    as.data.frame()
  # reorder rows to match manifest
  sub <- sub[match(manifest$gene, sub$gene), , drop = FALSE]
  out <- sub[, state_order_ext, drop = FALSE]
  colnames(out) <- state_short[state_order_ext]
  out
}

c1 <- build_chunk(samples_use[1], "mean_expr")       # SUR1090 expr
c2 <- build_chunk(samples_use[1], "pct_detected")    # SUR1090 pct
c3 <- build_chunk(samples_use[2], "mean_expr")       # SUR1072 expr
c4 <- build_chunk(samples_use[2], "pct_detected")    # SUR1072 pct

sep_col_df <- data.frame(Sep = rep("", nrow(manifest)), stringsAsFactors = FALSE)

final_df <- cbind(
  data.frame(State = as.character(manifest$panel), Gene = manifest$gene, stringsAsFactors = FALSE),
  c1, sep_col_df, c2, sep_col_df, c3, sep_col_df, c4
)

nc <- ncol(final_df)
nr <- nrow(final_df)

# Column index tracking
gene_col <- 2
c1_start <- 3;  c1_end <- 8    # 5 states + Total
sep1 <- 9
c2_start <- 10; c2_end <- 15
sep2 <- 16
c3_start <- 17; c3_end <- 22
sep3 <- 23
c4_start <- 24; c4_end <- 29

# --- Row 1: super headers ---
writeData(wb, sn, "State", startCol = 1, startRow = 1)
mergeCells(wb, sn, cols = 1, rows = 1:2)
writeData(wb, sn, "Gene", startCol = 2, startRow = 1)
mergeCells(wb, sn, cols = 2, rows = 1:2)

writeData(wb, sn, "SUR1090 Untreated — Expression", startCol = c1_start, startRow = 1)
mergeCells(wb, sn, cols = c1_start:c1_end, rows = 1)
addStyle(wb, sn, chunk1_hdr, rows = 1, cols = c1_start:c1_end, gridExpand = TRUE, stack = TRUE)

writeData(wb, sn, "SUR1090 Untreated — % Detected", startCol = c2_start, startRow = 1)
mergeCells(wb, sn, cols = c2_start:c2_end, rows = 1)
addStyle(wb, sn, chunk2_hdr, rows = 1, cols = c2_start:c2_end, gridExpand = TRUE, stack = TRUE)

writeData(wb, sn, "SUR1072 Untreated — Expression", startCol = c3_start, startRow = 1)
mergeCells(wb, sn, cols = c3_start:c3_end, rows = 1)
addStyle(wb, sn, chunk3_hdr, rows = 1, cols = c3_start:c3_end, gridExpand = TRUE, stack = TRUE)

writeData(wb, sn, "SUR1072 Untreated — % Detected", startCol = c4_start, startRow = 1)
mergeCells(wb, sn, cols = c4_start:c4_end, rows = 1)
addStyle(wb, sn, chunk4_hdr, rows = 1, cols = c4_start:c4_end, gridExpand = TRUE, stack = TRUE)

# --- Row 2: sub-headers (state short names) ---
row2 <- c("State", "Gene",
          state_short[state_order_ext], "",
          state_short[state_order_ext], "",
          state_short[state_order_ext], "",
          state_short[state_order_ext])
for (j in seq_along(row2)) writeData(wb, sn, row2[j], startCol = j, startRow = 2)
addStyle(wb, sn, header_style, rows = 2, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
addStyle(wb, sn, chunk1_hdr, rows = 2, cols = c1_start:c1_end, gridExpand = TRUE, stack = TRUE)
addStyle(wb, sn, chunk2_hdr, rows = 2, cols = c2_start:c2_end, gridExpand = TRUE, stack = TRUE)
addStyle(wb, sn, chunk3_hdr, rows = 2, cols = c3_start:c3_end, gridExpand = TRUE, stack = TRUE)
addStyle(wb, sn, chunk4_hdr, rows = 2, cols = c4_start:c4_end, gridExpand = TRUE, stack = TRUE)

# --- Data rows (starting row 3) ---
writeData(wb, sn, final_df, startCol = 1, startRow = 3, colNames = FALSE)

data_rows <- 3:(nr + 2)

# Gene bold
addStyle(wb, sn, gene_style, rows = data_rows, cols = gene_col, stack = TRUE)

# Expression formatting (3 decimals)
expr_cols <- c(c1_start:c1_end, c3_start:c3_end)
addStyle(wb, sn, num3, rows = data_rows, cols = expr_cols, gridExpand = TRUE, stack = TRUE)

# Pct formatting
pct_cols <- c(c2_start:c2_end, c4_start:c4_end)
addStyle(wb, sn, pct_fmt, rows = data_rows, cols = pct_cols, gridExpand = TRUE, stack = TRUE)

# Separator columns
for (sc in c(sep1, sep2, sep3)) {
  addStyle(wb, sn, sep_style, rows = 1:(nr + 2), cols = sc, stack = TRUE)
}

# Alternating rows
even_rows <- seq(4, nr + 2, by = 2)
if (length(even_rows) > 0) {
  addStyle(wb, sn, even_row_style, rows = even_rows, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
}

# State column: merge cells per state block and colour
sc_state_cols_full <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

state_runs <- rle(as.character(manifest$panel))
run_start <- 3
for (i in seq_along(state_runs$lengths)) {
  st <- state_runs$values[i]
  len <- state_runs$lengths[i]
  run_end <- run_start + len - 1
  mergeCells(wb, sn, cols = 1, rows = run_start:run_end)
  st_hex <- sc_state_cols_full[st]
  vert_style <- createStyle(
    textRotation = 90, halign = "center", valign = "center",
    textDecoration = "bold", fontColour = "#FFFFFF", fgFill = st_hex,
    border = "TopBottomLeftRight"
  )
  addStyle(wb, sn, vert_style, rows = run_start:run_end, cols = 1, stack = TRUE)
  # Gene text colour
  gene_col_style <- createStyle(textDecoration = "bold", fontName = "Consolas", fontColour = st_hex)
  addStyle(wb, sn, gene_col_style, rows = run_start:run_end, cols = 2, stack = TRUE)
  run_start <- run_end + 1
}

# Conditional formatting — expression: white -> pink -> dark red
if (nr > 1) {
  for (cc in expr_cols) {
    conditionalFormatting(wb, sn, cols = cc, rows = data_rows,
      style = c("#FFFFFF", "#FB8A8A", "#B22222"),
      rule = c(0, 1.5, 3.5), type = "colourScale")
  }
  for (cc in pct_cols) {
    conditionalFormatting(wb, sn, cols = cc, rows = data_rows,
      style = c("#FFFFFF", "#27AE60"),
      rule = c(0, 1), type = "colourScale")
  }
}

# Borders at chunk boundaries
for (bc in c(gene_col, c1_start, sep1, c2_start, sep2, c3_start, sep3, c4_start)) {
  addStyle(wb, sn, border_left, rows = 1:(nr + 2), cols = bc, stack = TRUE)
}

# Column widths
setColWidths(wb, sn, cols = 1, widths = 5)
setColWidths(wb, sn, cols = 2, widths = 14)
setColWidths(wb, sn, cols = c(c1_start:c1_end, c2_start:c2_end, c3_start:c3_end, c4_start:c4_end), widths = 11)
setColWidths(wb, sn, cols = c(sep1, sep2, sep3), widths = 2)

freezePane(wb, sn, firstActiveRow = 3, firstActiveCol = 3)

out_xlsx <- file.path(out_dir, "Auto_marker_sample_expression.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved Excel to: ", out_xlsx)

####################
# ============ DOT PLOT PDF ============
####################
message("Building specificity dot plot...")

plot_df <- res_df %>%
  mutate(
    state = factor(state, levels = state_order_ext),
    sample_short = sample_labels_short[sample],
    sample_short = factor(sample_short, levels = c("SUR1090 Untreated", "SUR1072 Untreated")),
    panel = state  # each gene belongs to its own state panel
  )

# Gene order: reversed manifest order (bottom to top)
gene_order <- rev(manifest$gene)
plot_df$gene <- factor(plot_df$gene, levels = gene_order)

# Panel membership for faceting
plot_df <- plot_df %>%
  left_join(manifest %>% rename(facet_panel = panel), by = "gene") %>%
  mutate(facet_panel = factor(facet_panel, levels = state_order))

# Normalise expression within each sample per gene (0-1 scaling)
# The range (0 and 1) is defined only by the 5 state values; Total is then mapped onto this scale.
plot_df <- plot_df %>%
  group_by(gene, sample) %>%
  mutate(
    norm_expr = {
      state_vals <- mean_expr[state != "Total"]
      rng <- range(state_vals, na.rm = TRUE)
      if (all(is.na(state_vals)) || diff(rng) == 0) {
        rep(0.5, n())
      } else {
        (mean_expr - rng[1]) / diff(rng)
      }
    }
  ) %>%
  ungroup()

# Build x-axis as state × sample interaction
plot_df$x_label <- paste0(as.character(plot_df$state), "\n", as.character(plot_df$sample_short))

# For proper ordering within facets
plot_df$x_pos <- interaction(plot_df$state, plot_df$sample_short, lex.order = TRUE)

det_plot <- ggplot(plot_df, aes(x = sample_short, y = gene)) +
  geom_point(aes(size = pct_detected, fill = norm_expr), shape = 21, color = "grey35", stroke = 0.35) +
  facet_grid(
    facet_panel ~ state,
    scales = "free",
    space = "free",
    switch = "y",
    labeller = labeller(
      facet_panel = setNames(
        sprintf("<span style='color:%s'>%s</span>", state_cols[state_order], state_order),
        state_order
      ),
      state = setNames(
        c(sprintf("<span style='color:%s'>%s</span>", state_cols[state_order], state_order),
          "<b>Total</b>"),
        state_order_ext
      )
    )
  ) +
  scale_size_continuous(range = c(0.5, 5.5), labels = scales::percent_format(accuracy = 1), name = "% Detected") +
  scale_fill_gradient2(
    low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0.5,
    name = "Normalised\nExpression",
    limits = c(0, 1)
  ) +
  scale_x_discrete(labels = function(x) gsub(" Untreated", "\nUntreated", x)) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(size = 7, angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8),
    strip.placement = "outside",
    strip.text.y.left = ggtext::element_markdown(angle = 0, face = "bold", size = 7.5),
    strip.text.x = ggtext::element_markdown(face = "bold", size = 8, lineheight = 1.1),
    panel.grid.minor = element_blank(),
    panel.spacing.x = unit(2.5, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  ) +
  labs(x = NULL, y = NULL)

# Save the plot
out_pdf <- file.path(out_dir, "Auto_marker_sample_specificity_dotplot.pdf")

ggsave(out_pdf, det_plot, width = 14, height = 6.5, useDingbats = FALSE)
message("Saved PDF to: ", out_pdf)

message("Done.")
