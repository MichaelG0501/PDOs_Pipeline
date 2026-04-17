####################
# Auto_PDO_scAtlas_scenic_comparison.R
# Compare SCENIC regulon activities (RSS) between scAtlas and PDOs.
####################
# NOTE ON SCORES:
# The prompt asked to confirm if the score is comparable between datasets (like AUCell).
# AUCell scores evaluate the recovery of gene targets within the expression ranking of each cell,
# yielding a metric that is generally comparable across datasets.
# However, to highlight *specificity* to cell states (which regulons mark which state), the pipeline 
# computes the Regulon Specificity Score (RSS) based on the Jensen-Shannon divergence of AUCs. 
# We use RSS here to capture the specific regulatory profile per MP/State as it is better suited for 
# distinguishing states. To ensure both datasets visually align on the same relative scale and emphasize 
# structural similarities, the combined RSS matrix is row-scaled (Z-scored) before hierarchical clustering.
####################

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

sc_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/final_mp_scenic"
pdo_dir <- "final_mp_scenic"

# 1. Load RSS matrices
sc_mp_rss <- readRDS(file.path(sc_dir, "Auto_final_mp_scenic_rss.rds"))
sc_st_rss <- readRDS(file.path(sc_dir, "Auto_final_mp_scenic_state_rss.rds"))

pdo_mp_rss <- readRDS(file.path(pdo_dir, "Auto_PDO_final_mp_scenic_rss.rds"))
pdo_st_rss <- readRDS(file.path(pdo_dir, "Auto_PDO_final_mp_scenic_state_rss.rds"))

# 2. Load raw AUC matrices and cell metadata for state-level AUC averages
message("Loading cell-level AUC matrices and metadata...")
sc_auc_mat <- readRDS(file.path(sc_dir, "Auto_final_mp_scenic_regulon_auc.rds"))
pdo_auc_mat <- readRDS(file.path(pdo_dir, "Auto_PDO_final_mp_scenic_regulon_auc.rds"))

sc_metadata <- read.csv(file.path(sc_dir, "Auto_final_mp_scenic_selected_cells.csv"), stringsAsFactors=FALSE)
pdo_metadata <- read.csv(file.path(pdo_dir, "Auto_PDO_final_mp_scenic_selected_cells.csv"), stringsAsFactors=FALSE)

# 3. Load cell mappings to get MP_group colours
sc_selected <- read.csv(file.path(sc_dir, "Auto_final_mp_scenic_selected_cells.csv"), stringsAsFactors=FALSE)
pdo_selected <- read.csv(file.path(pdo_dir, "Auto_PDO_final_mp_scenic_selected_cells.csv"), stringsAsFactors=FALSE)

sc_mp_anno <- sc_selected %>% distinct(final_mp_label, mp_group) %>% mutate(Dataset = "scAtlas")
pdo_mp_anno <- pdo_selected %>% distinct(final_mp_label, mp_group) %>% mutate(Dataset = "PDO")

format_regulon_name <- function(x) {
  x <- gsub(" \\([0-9]+g\\)$", "", x)
  x <- gsub(" \\([0-9]+ genes\\)$", "", x)
  gsub("_extended$", "", x)
}

clean_mat <- function(mat) {
  rownames(mat) <- format_regulon_name(rownames(mat))
  # Average duplicated rows if any extend/normal motifs collapse to same TF
  if (any(duplicated(rownames(mat)))) {
    rsum <- rowsum(mat, rownames(mat))
    rcount <- table(rownames(mat))
    mat <- rsum / as.numeric(rcount[rownames(rsum)])
  }
  mat
}

# Fix rownames to base regulon names to allow precise intersection
sc_mp_rss <- clean_mat(sc_mp_rss)
pdo_mp_rss <- clean_mat(pdo_mp_rss)
sc_st_rss <- clean_mat(sc_st_rss)
pdo_st_rss <- clean_mat(pdo_st_rss)

####################
# Specificity Gap Logic
####################
# Replace RSS with Specificity Gap: RSS(state) - max(RSS(other states))
calc_rss_gap <- function(rss_mat, states_of_interest) {
  # Subset to defined states for fair comparison
  mat <- rss_mat[, intersect(states_of_interest, colnames(rss_mat)), drop=FALSE]
  gap_mat <- matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
  rownames(gap_mat) <- rownames(mat)
  colnames(gap_mat) <- colnames(mat)
  
  for (i in 1:nrow(mat)) {
    row_vals <- mat[i, ]
    for (j in 1:ncol(mat)) {
      gap_mat[i, j] <- row_vals[j] - max(row_vals[-j], na.rm=TRUE)
    }
  }
  gap_mat
}

# Define states of interest for gap calculation
sc_defined_states <- c("Classic Proliferative", "Basal to Intestinal Metaplasia", "Stress-adaptive", "SMG-like Metaplasia", "Immune Infiltrating", "3CA_EMT_and_Protein_maturation")
pdo_defined_states <- c("Classic Proliferative", "Basal to Intest. Meta", "Stress-adaptive", "SMG-like Metaplasia", "3CA_EMT_and_Protein_maturation")

sc_st_gap  <- calc_rss_gap(sc_st_rss, sc_defined_states)
pdo_st_gap <- calc_rss_gap(pdo_st_rss, pdo_defined_states)

# Update Excel helper to use Gap
get_state_gap_vec <- function(gap_mat, state_name, all_regs) {
  if (!is.null(gap_mat) && state_name %in% colnames(gap_mat)) {
    rn <- format_regulon_name(rownames(gap_mat))
    df <- data.frame(Regulon = rn, Gap = gap_mat[, state_name], stringsAsFactors = FALSE)
    df <- df %>% group_by(Regulon) %>% summarize(Gap = mean(Gap, na.rm=TRUE), .groups="drop")
    return(df$Gap[match(all_regs, df$Regulon)])
  }
  return(rep(NA_real_, length(all_regs)))
}

run_heatmap <- function(mat1, mat2, meta1, meta2, out_pdf, title, type="MP") {
  com_regs <- intersect(rownames(mat1), rownames(mat2))
  message(sprintf("Found %d common regulons for %s", length(com_regs), type))

  # Combine matrices using intersected regulons
  mat_comb <- cbind(mat1[com_regs, , drop=FALSE], mat2[com_regs, , drop=FALSE])

  # Row-scale to emphasize relative differences strictly across the integrated set
  mat_scaled <- t(scale(t(mat_comb)))
  mat_scaled[!is.finite(mat_scaled)] <- 0

  # Build annotations
  if (type == "MP") {
    meta_comb <- bind_rows(meta1, meta2)
    meta_comb <- meta_comb[match(colnames(mat_comb), meta_comb$final_mp_label), ]
    
    group_cols <- c(
      "Cell cycle" = "#D4AF37",
      "Cell Cycle" = "#D4AF37",
      "Classic Proliferative" = "#E41A1C",
      "Basal to Intest. Meta" = "#4DAF4A",
      "Basal to Intestinal Metaplasia" = "#4DAF4A",
      "Stress-adaptive" = "#984EA3",
      "SMG-like Metaplasia" = "#FF7F00",
      "Pan-cancer 3CA" = "#6A3D9A",
      "Immune Infiltrating" = "#A65628",
      "Other" = "grey70"
    )
    ha <- HeatmapAnnotation(
      Dataset = meta_comb$Dataset,
      Group = meta_comb$mp_group,
      col = list(
        Dataset = c("scAtlas" = "grey30", "PDO" = "grey80"),
        Group = group_cols
      ),
      show_annotation_name = TRUE
    )
  } else {
    meta_comb <- data.frame(
      State = colnames(mat_comb),
      Dataset = c(rep("scAtlas", ncol(mat1)), rep("PDO", ncol(mat2)))
    )

    state_cols <- c(
      "Classic Proliferative" = "#E41A1C",
      "Basal to Intest. Meta" = "#4DAF4A",
      "Basal to Intestinal Metaplasia" = "#4DAF4A",
      "Stress-adaptive" = "#984EA3",
      "SMG-like Metaplasia" = "#FF7F00",
      "3CA_EMT_and_Protein_maturation" = "#377EB8",
      "Immune Infiltrating" = "#A65628",
      "Unresolved" = "grey80",
      "Hybrid" = "black"
    )
    
    ha <- HeatmapAnnotation(
      Dataset = meta_comb$Dataset,
      State = meta_comb$State,
      col = list(
        Dataset = c("scAtlas" = "grey30", "PDO" = "grey80"),
        State = state_cols
      ),
      show_annotation_name = TRUE
    )
  }
  
  # Sequential professional red scale for Z-scored RSS
  col_fun <- colorRamp2(c(0, 1.25, 2.5), c("#FFFFFF", "#FB8A8A", "#B22222"))
  
  pdf(out_pdf, width = 18, height = 15, useDingbats = FALSE)
  draw(
    Heatmap(
      mat_scaled,
      name = "Scaled RSS",
      col = col_fun,
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_column_dend = TRUE,
      row_names_side = "left",
      row_names_gp = gpar(fontsize = max(4, min(8, 600/length(com_regs)))),
      column_names_gp = gpar(fontsize = 10),
      column_names_rot = 45,
      show_row_names = TRUE
    ),
    merge_legend = TRUE,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )
  grid.text(
    title,
    x = unit(4, "mm"),
    y = unit(1, "npc") - unit(4, "mm"),
    just = c("left", "top"),
    gp = gpar(fontsize = 15, fontface = "bold")
  )
  dev.off()
}

out_dir_plot <- file.path(pdo_dir, "scAtlas_comparison")
dir.create(out_dir_plot, showWarnings = FALSE)

# Heatmap for MPs
run_heatmap(
  sc_mp_rss, pdo_mp_rss, 
  sc_mp_anno, pdo_mp_anno, 
  file.path(out_dir_plot, "Auto_scenic_comparison_MP_heatmap.pdf"),
  "SCENIC Regulon Specificity Comparison (scAtlas vs PDO MPs)",
  type="MP"
)

# Heatmap for States
run_heatmap(
  sc_st_rss, pdo_st_rss, 
  NULL, NULL, 
  file.path(out_dir_plot, "Auto_scenic_comparison_State_heatmap.pdf"),
  "SCENIC Regulon Specificity Comparison (scAtlas vs PDO States)",
  type="State"
)
message("Saved comparison plots to ", out_dir_plot)

####################
# Generate Comparative Excel Table
####################
library(openxlsx)

message("Generating comparative Excel table...")

# 1. Start building Overview based on alphabetical unique regulon names
raw_regs <- unique(c(rownames(sc_st_rss), rownames(pdo_st_rss)))
all_regulons <- sort(unique(format_regulon_name(raw_regs)))
overview_df <- data.frame(Regulon = all_regulons, stringsAsFactors = FALSE)

# Shared state mapping
shared_states <- list(
  "Classic Proliferative" = list(sc = "Classic Proliferative", pdo = "Classic Proliferative"),
  "Stress-adaptive" = list(sc = "Stress-adaptive", pdo = "Stress-adaptive"),
  "Basal Metaplasia" = list(sc = "Basal to Intestinal Metaplasia", pdo = "Basal to Intest. Meta"),
  "SMG-like Metaplasia" = list(sc = "SMG-like Metaplasia", pdo = "SMG-like Metaplasia"),
  "3CA EMT" = list(sc = "3CA_EMT_and_Protein_maturation", pdo = "3CA_EMT_and_Protein_maturation")
)

# Helper functions for state-level aggregation
get_state_rss <- function(rss_mat, state_name) {
  if (!is.null(rss_mat) && state_name %in% colnames(rss_mat)) {
    rn <- format_regulon_name(rownames(rss_mat))
    df <- data.frame(Regulon = rn, RSS = rss_mat[, state_name], stringsAsFactors = FALSE)
    return(df %>% group_by(Regulon) %>% summarize(RSS = mean(RSS, na.rm=TRUE), .groups="drop"))
  }
  return(data.frame(Regulon = character(), RSS = numeric(), stringsAsFactors = FALSE))
}

get_state_auc <- function(auc_obj, metadata, state_name) {
  auc_mat <- AUCell::getAUC(auc_obj)
  cells <- metadata$cell[metadata$final_state == state_name]
  cells <- intersect(cells, colnames(auc_mat))
  if (length(cells) > 0) {
    rn <- format_regulon_name(rownames(auc_mat))
    means <- rowMeans(as.matrix(auc_mat[, cells, drop=FALSE]))
    df <- data.frame(Regulon = rn, AUC = means, stringsAsFactors = FALSE)
    return(df %>% group_by(Regulon) %>% summarize(AUC = mean(AUC, na.rm=TRUE), .groups="drop"))
  }
  return(data.frame(Regulon = character(), AUC = numeric(), stringsAsFactors = FALSE))
}

## 1. Initialize result objects
raw_regs <- unique(c(rownames(sc_st_rss), rownames(pdo_st_rss)))
all_regulons <- sort(unique(format_regulon_name(raw_regs)))
overview_df <- data.frame(Regulon = all_regulons, stringsAsFactors = FALSE)

state_header_pos <- list()
current_col <- 2

# Shared state mapping
shared_states <- list(
  "Classic Proliferative" = list(sc = "Classic Proliferative", pdo = "Classic Proliferative"),
  "Stress-adaptive" = list(sc = "Stress-adaptive", pdo = "Stress-adaptive"),
  "Basal Metaplasia" = list(sc = "Basal to Intestinal Metaplasia", pdo = "Basal to Intest. Meta"),
  "SMG-like Metaplasia" = list(sc = "SMG-like Metaplasia", pdo = "SMG-like Metaplasia"),
  "3CA EMT" = list(sc = "3CA_EMT_and_Protein_maturation", pdo = "3CA_EMT_and_Protein_maturation")
)

# Helper to get specific state data
get_state_rss_vec <- function(rss_mat, state_name) {
  df <- get_state_rss(rss_mat, state_name)
  res <- df$RSS[match(all_regulons, df$Regulon)]
  return(res)
}

# --- Section A: Interleaved RSS (Original Style) ---
for (st_label in names(shared_states)) {
  sc_st <- shared_states[[st_label]]$sc
  pdo_st <- shared_states[[st_label]]$pdo
  
  col_sc <- paste0(st_label, "_scAtlas")
  col_pdo <- paste0(st_label, "_PDO")
  col_comb <- paste0(st_label, "_Combined")
  
  overview_df[[col_sc]] <- get_state_gap_vec(sc_st_gap, sc_st, all_regulons)
  overview_df[[col_pdo]] <- get_state_gap_vec(pdo_st_gap, pdo_st, all_regulons)
  overview_df[[col_comb]] <- rowMeans(cbind(overview_df[[col_sc]], overview_df[[col_pdo]]), na.rm=TRUE)
  
  state_header_pos[[st_label]] <- list(start = current_col, middle = current_col + 1, end = current_col + 2)
  current_col <- current_col + 3
}

if ("Immune Infiltrating" %in% colnames(sc_st_rss)) {
  st_label <- "Immune Infiltrating"
  col_sc <- paste0(st_label, "_scAtlas")
  overview_df[[col_sc]] <- get_state_gap_vec(sc_st_gap, "Immune Infiltrating", all_regulons)
  state_header_pos[[st_label]] <- list(start = current_col, middle = current_col, end = current_col)
  current_col <- current_col + 1
}

# --- Section B: Separator ---
overview_df$Sep <- ""
sep_col_idx <- current_col
current_col <- current_col + 1

# --- Section C: Grouped AUC (Marker style) ---
state_abbrev <- c(
  "Classic Proliferative" = "ClassProlif",
  "Stress-adaptive" = "StressAdapt",
  "Basal Metaplasia" = "BasalMeta",
  "SMG-like Metaplasia" = "SMG-like",
  "3CA EMT" = "3CA_EMT",
  "Immune Infiltrating" = "ImmuneInfil"
)

get_state_auc_vec <- function(auc_obj, metadata, state_name) {
  df <- get_state_auc(auc_obj, metadata, state_name)
  res <- df$AUC[match(all_regulons, df$Regulon)]
  return(res)
}

auc_sc_start <- current_col
# scAtlas AUC Group
for (st in names(state_abbrev)) {
  # Map shared state to sc name, or use name directly for immune
  sc_name <- if (st %in% names(shared_states)) shared_states[[st]]$sc else st
  if (sc_name %in% sc_metadata$final_state) {
    col_name <- paste0("sc_AUC_", state_abbrev[st])
    overview_df[[col_name]] <- get_state_auc_vec(sc_auc_mat, sc_metadata, sc_name)
    current_col <- current_col + 1
  }
}
auc_sc_end <- current_col - 1

auc_pdo_start <- current_col
# PDO AUC Group
for (st in names(state_abbrev)) {
  if (st %in% names(shared_states)) {
    pdo_name <- shared_states[[st]]$pdo
    if (pdo_name %in% pdo_metadata$final_state) {
       col_name <- paste0("pdo_AUC_", state_abbrev[st])
       overview_df[[col_name]] <- get_state_auc_vec(pdo_auc_mat, pdo_metadata, pdo_name)
       current_col <- current_col + 1
    }
  }
}
auc_pdo_end <- current_col - 1

# Values for color scale
# Gap scale - use fixed professional diverging range
min_gap <- -0.2
max_gap <- 0.2

auc_cols_names <- colnames(overview_df)[(sep_col_idx+1):ncol(overview_df)]
valid_auc_vals <- as.matrix(overview_df[, auc_cols_names])[is.finite(as.matrix(overview_df[, auc_cols_names]))]
min_auc <- if(length(valid_auc_vals) > 0) quantile(valid_auc_vals, 0.05, na.rm=TRUE) else 0
max_auc <- if(length(valid_auc_vals) > 0) quantile(valid_auc_vals, 0.95, na.rm=TRUE) else 0.5

# Display names for Row 2
display_colnames <- colnames(overview_df)
# For RSS, show dataset names (scAtlas, PDO, Combined)
display_colnames[grepl("_scAtlas$", display_colnames)] <- "scAtlas"
display_colnames[grepl("_PDO$", display_colnames)] <- "PDO"
display_colnames[grepl("_Combined$", display_colnames)] <- "Combined"

display_colnames[colnames(overview_df) == "Sep"] <- ""

# For AUC, keep abbreviations
display_colnames[grepl("^sc_AUC_|^pdo_AUC_", colnames(overview_df))] <- gsub("^sc_AUC_|^pdo_AUC_", "", colnames(overview_df)[grepl("^sc_AUC_|^pdo_AUC_", colnames(overview_df))])

# 2. State-Specific Sheets (Original 3-col style)
state_sheets <- list()
for (st_label in names(shared_states)) {
  sc_col <- paste0(st_label, "_scAtlas")
  pdo_col <- paste0(st_label, "_PDO")
  comb_col <- paste0(st_label, "_Combined")
  st_df <- overview_df[, c("Regulon", sc_col, pdo_col, comb_col)]
  colnames(st_df) <- c("Regulon", "scAtlas", "PDO", "Combined")
  state_sheets[[st_label]] <- st_df %>% arrange(desc(Combined))
}

# 3. Write Excel
wb <- createWorkbook()

# Match styles from Auto_marker_comparison_excel.R
sc_header_style <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  fontColour = "#FFFFFF", fgFill = "#2C3E50", border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
pdo_header_style <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  fontColour = "#FFFFFF", fgFill = "#8E44AD", border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
comb_header_style <- createStyle(
  textDecoration = "bold", halign = "center", valign = "center",
  fontColour = "#FFFFFF", fgFill = "#27AE60", border = "Bottom", borderStyle = "medium", wrapText = TRUE
)
boldStyle <- createStyle(textDecoration = "bold", halign = "center", valign = "center", border = "Bottom", borderStyle = "thick")
numStyle  <- createStyle(numFmt = "0.000")
sep_style <- createStyle(fgFill = "#D5D8DC", border = "LeftRight", borderColour = "#95A5A6")
gene_style <- createStyle(textDecoration = "bold", fontName = "Consolas")
border_style <- createStyle(border = "Left", borderStyle = "medium", borderColour = "#2C3E50")

# --- Combined Overview Sheet ---
addWorksheet(wb, "Combined Overview")
writeData(wb, "Combined Overview", "Regulon", startCol=1, startRow=1)
mergeCells(wb, "Combined Overview", cols = 1, rows = 1:2)

# Row 1 RSS Labels (Interleaved style)
for (st_label in names(state_header_pos)) {
  pos <- state_header_pos[[st_label]]
  writeData(wb, "Combined Overview", st_label, startCol=pos$middle, startRow=1)
  addStyle(wb, "Combined Overview", border_style, rows = 1:(nrow(overview_df)+2), cols = pos$start, stack = TRUE)
}

# Row 1 AUC Labels (Grouped style)
writeData(wb, "Combined Overview", "scAtlas AUC", startCol = auc_sc_start, startRow = 1)
writeData(wb, "Combined Overview", "PDO AUC",     startCol = auc_pdo_start, startRow = 1)
addStyle(wb, "Combined Overview", sc_header_style, rows = 1, cols = auc_sc_start:auc_sc_end, gridExpand = TRUE)
addStyle(wb, "Combined Overview", pdo_header_style, rows = 1, cols = auc_pdo_start:auc_pdo_end, gridExpand = TRUE)

# Row 2 Labels
writeData(wb, "Combined Overview", t(display_colnames[-1]), startCol=2, startRow=2, colNames=FALSE)
addStyle(wb, "Combined Overview", boldStyle, rows=1:2, cols=1:ncol(overview_df), gridExpand=TRUE, stack=TRUE)

# Data
writeData(wb, "Combined Overview", overview_df, startCol=1, startRow=3, colNames=FALSE)
addStyle(wb, "Combined Overview", gene_style, rows = 3:(nrow(overview_df)+2), cols = 1, stack = TRUE)
addStyle(wb, "Combined Overview", numStyle, rows = 3:(nrow(overview_df)+2), cols = 2:ncol(overview_df), gridExpand = TRUE, stack=TRUE)
addStyle(wb, "Combined Overview", sep_style, rows = 1:(nrow(overview_df)+2), cols = sep_col_idx, stack = TRUE)

conditionalFormatting(wb, "Combined Overview", cols = 2:(sep_col_idx-1), rows = 3:(nrow(overview_df)+2), 
                      style = c("#1D4E89", "#F8F4EC", "#B22222"), rule = c(min_gap, 0, max_gap), type = "colourScale")

# Color scale - AUC section (Requested 5th to 95th)
conditionalFormatting(wb, "Combined Overview", cols = (sep_col_idx+1):ncol(overview_df), rows = 3:(nrow(overview_df)+2), 
                      style = c("#FFFFFF", "#FB8A8A", "#B22222"), rule = c(min_auc, (min_auc + max_auc)/2, max_auc), type = "colourScale")

# Left border at AUC section boundaries
addStyle(wb, "Combined Overview", border_style, rows = 1:(nrow(overview_df)+2), cols = auc_sc_start, stack = TRUE)
addStyle(wb, "Combined Overview", border_style, rows = 1:(nrow(overview_df)+2), cols = auc_pdo_start, stack = TRUE)
addStyle(wb, "Combined Overview", border_style, rows = 1:(nrow(overview_df)+2), cols = sep_col_idx, stack = TRUE)

# Column widths
setColWidths(wb, "Combined Overview", cols = 1, widths = 25)
setColWidths(wb, "Combined Overview", cols = 2:ncol(overview_df), widths = 12)
setColWidths(wb, "Combined Overview", cols = sep_col_idx, widths = 4)
freezePane(wb, "Combined Overview", firstActiveRow = 3, firstActiveCol = 2)

# --- State Sheets ---
for (st_label in names(state_sheets)) {
  sheet_name <- substr(st_label, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, state_sheets[[st_label]], headerStyle = boldStyle)
  addStyle(wb, sheet_name, numStyle, rows=2:(nrow(state_sheets[[st_label]])+1), cols=2:4, gridExpand=TRUE, stack=TRUE)
  # Separate scale for per-sheet data
  sheet_vals <- as.matrix(state_sheets[[st_label]][, 2:4])
  sheet_max <- quantile(sheet_vals[is.finite(sheet_vals)], 0.95, na.rm=TRUE)
  if (is.na(sheet_max) || sheet_max == 0) sheet_max <- 0.5
  
  conditionalFormatting(wb, sheet_name, cols = 2:4, rows = 2:(nrow(state_sheets[[st_label]])+1), 
                        style = c("#1D4E89", "#F8F4EC", "#B22222"), rule = c(-0.2, 0, 0.2), type = "colourScale")
  setColWidths(wb, sheet_name, cols = 1, widths = 25)
  setColWidths(wb, sheet_name, cols = 2:4, widths = 12)
  freezePane(wb, sheet_name, firstActiveRow = 2, firstActiveCol = 2)
}

out_xlsx <- file.path(pdo_dir, "Auto_scRef_PDO_scenic_comparison.xlsx")
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Saved beautiful comparative Excel table to ", out_xlsx)

# Generate 3-page RSS Heatmap PDF
####################
message("Generating 3-page RSS heatmap PDF...")

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

sc_state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "Immune Infiltrating" = "#377EB8",
  "3CA_EMT_and_Protein_maturation" = "#666666"
)

pdo_state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive" = "#984EA3",
  "SMG-like Metaplasia" = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8"
)

# Color scale for Specificity Gap
col_fun <- colorRamp2(c(-0.25, 0, 0.25), c("#1D4E89", "#F8F4EC", "#B22222"))

# Heatmap 1: scAtlas 6 states
sc_hm_list <- list()
for(st in sc_state_order) {
  top5 <- names(sort(sc_st_gap[, st], decreasing=TRUE)[1:5])
  sc_hm_list[[st]] <- data.frame(Regulon = top5, State = st, stringsAsFactors=FALSE)
}
sc_hm_df <- do.call(rbind, sc_hm_list)
sc_hm_df$GapVal <- sc_st_gap[cbind(sc_hm_df$Regulon, sc_hm_df$State)]
sc_hm_df <- sc_hm_df[order(sc_hm_df$GapVal, decreasing=TRUE), ]
sc_hm_df <- sc_hm_df[!duplicated(sc_hm_df$Regulon), ]
sc_hm_df$State <- factor(sc_hm_df$State, levels = sc_state_order)
sc_hm_df <- sc_hm_df[order(sc_hm_df$State), ]

sc_hm_regs <- sc_hm_df$Regulon
sc_reg_state <- sc_hm_df$State

sc_hm_plot <- sc_st_gap[sc_hm_regs, sc_state_order, drop=FALSE]

row_ann_sc <- rowAnnotation(
  State = sc_reg_state,
  col = list(State = sc_state_cols),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

top_ann_sc <- HeatmapAnnotation(
  Top_State = factor(sc_state_order, levels = sc_state_order),
  col = list(Top_State = sc_state_cols),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(4, "mm")
)

ht_sc <- Heatmap(
  sc_hm_plot, name = "Specificity\nGap",
  top_annotation = top_ann_sc, left_annotation = row_ann_sc, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = sc_reg_state, row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, border = TRUE
)


# Heatmap 2: PDO 5 states
pdo_hm_list <- list()
for(st in pdo_state_order) {
  top5 <- names(sort(pdo_st_gap[, st], decreasing=TRUE)[1:5])
  pdo_hm_list[[st]] <- data.frame(Regulon = top5, State = st, stringsAsFactors=FALSE)
}
pdo_hm_df <- do.call(rbind, pdo_hm_list)
pdo_hm_df$GapVal <- pdo_st_gap[cbind(pdo_hm_df$Regulon, pdo_hm_df$State)]
pdo_hm_df <- pdo_hm_df[order(pdo_hm_df$GapVal, decreasing=TRUE), ]
pdo_hm_df <- pdo_hm_df[!duplicated(pdo_hm_df$Regulon), ]
pdo_hm_df$State <- factor(pdo_hm_df$State, levels = pdo_state_order)
pdo_hm_df <- pdo_hm_df[order(pdo_hm_df$State), ]

pdo_hm_regs <- pdo_hm_df$Regulon
pdo_reg_state <- pdo_hm_df$State

pdo_hm_plot <- pdo_st_gap[pdo_hm_regs, pdo_state_order, drop=FALSE]

row_ann_pdo <- rowAnnotation(
  State = pdo_reg_state,
  col = list(State = pdo_state_cols),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

top_ann_pdo <- HeatmapAnnotation(
  Top_State = factor(pdo_state_order, levels = pdo_state_order),
  col = list(Top_State = pdo_state_cols),
  show_annotation_name = FALSE, show_legend = FALSE,
  simple_anno_size = unit(4, "mm")
)

ht_pdo <- Heatmap(
  pdo_hm_plot, name = "Specificity\nGap",
  top_annotation = top_ann_pdo, left_annotation = row_ann_pdo, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = pdo_reg_state, row_title_rot = 0,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, border = TRUE
)


# Heatmap 3: Combined Heatmap
comb_regs_list <- list()
for (st_label in names(shared_states)) {
  sc_st <- shared_states[[st_label]]$sc
  pdo_st <- shared_states[[st_label]]$pdo
  
  common <- intersect(rownames(sc_st_gap), rownames(pdo_st_gap))
  sc_vals <- sc_st_gap[common, sc_st]
  pdo_vals <- pdo_st_gap[common, pdo_st]
  
  has_support <- sc_vals > 0 & pdo_vals > 0
  if (any(has_support)) {
    comb_score <- (sc_vals[has_support] + pdo_vals[has_support]) / 2
    top5 <- names(sort(comb_score, decreasing=TRUE)[1:min(5, length(comb_score))])
    comb_regs_list[[sc_st]] <- data.frame(Regulon = top5, State = sc_st, Score = comb_score[top5], stringsAsFactors=FALSE)
  }
}

# Handle Immune Infiltrating separately (only scAtlas Gap)
if ("Immune Infiltrating" %in% sc_state_order) {
  im_regs <- names(sort(sc_st_gap[, "Immune Infiltrating"], decreasing=TRUE)[1:5])
  comb_regs_list[["Immune Infiltrating"]] <- data.frame(Regulon = im_regs, State = "Immune Infiltrating", Score = sc_st_gap[im_regs, "Immune Infiltrating"], stringsAsFactors=FALSE)
}

comb_regs_df <- do.call(rbind, comb_regs_list)
comb_regs_df <- comb_regs_df[order(comb_regs_df$Score, decreasing=TRUE), ]
comb_regs_df <- comb_regs_df[!duplicated(comb_regs_df$Regulon), ]
comb_regs_df$State <- factor(comb_regs_df$State, levels = sc_state_order)
comb_regs_df <- comb_regs_df[order(comb_regs_df$State), ]
comb_regs <- comb_regs_df$Regulon

mat_sc <- sc_st_gap[comb_regs, sc_state_order, drop=FALSE]
colnames(mat_sc) <- paste0("scAtlas::", sc_state_order)

mat_pdo <- matrix(NA, nrow=length(comb_regs), ncol=length(pdo_state_order))
rownames(mat_pdo) <- comb_regs
colnames(mat_pdo) <- paste0("PDO::", pdo_state_order)
regs_in_pdo <- comb_regs %in% rownames(pdo_st_gap)
mat_pdo[regs_in_pdo, ] <- pdo_st_gap[comb_regs[regs_in_pdo], pdo_state_order, drop=FALSE]

comb_plot <- cbind(mat_sc, mat_pdo)

row_ann_comb <- rowAnnotation(
  State = comb_regs_df$State,
  col = list(State = sc_state_cols),
  show_annotation_name = FALSE,
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
  comb_plot, name = "Specificity\nGap",
  top_annotation = top_ann_comb, left_annotation = row_ann_comb, col = col_fun,
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = comb_regs_df$State, row_title_rot = 0,
  column_split = col_split_factors,
  row_names_gp = gpar(fontsize = 8, fontface = "bold"),
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_names_rot = 45, border = TRUE
)

pdf(file.path(out_dir_plot, "Auto_scenic_comparison_RSS_heatmaps.pdf"), width = 17, height = 12, useDingbats = FALSE)
draw(ht_sc)
draw(ht_pdo)
draw(ht_comb)
dev.off()

####################
# Output Separate Top 5 Excel
####################
# Build a separate top5 Excel with 3 sheets reusing the main overview formatting
wb2 <- createWorkbook()

build_top5_sheet <- function(wb, sheet_name, mapping_df) {
  addWorksheet(wb, sheet_name)
  
  # Ensure ordered by state
  mapping_df <- mapping_df[order(mapping_df$State), ]
  regs <- mapping_df$Regulon
  
  # Prep sub_df from overview
  sub_df_base <- overview_df[overview_df$Regulon %in% regs, ]
  sub_df_base <- sub_df_base[match(regs, sub_df_base$Regulon), ]
  
  # Build final data with empty rows
  final_df <- data.frame()
  row_colors <- character()
  is_data_row <- logical()
  
  levels_present <- levels(mapping_df$State)
  for (i in seq_along(levels_present)) {
    st <- levels_present[i]
    st_regs <- mapping_df$Regulon[mapping_df$State == st]
    if (length(st_regs) == 0) next
    
    st_data <- sub_df_base[sub_df_base$Regulon %in% st_regs, ]
    final_df <- rbind(final_df, st_data)
    
    # Get color
    st_col <- sc_state_cols[st]
    if (is.na(st_col)) st_col <- "#000000"
    
    row_colors <- c(row_colors, rep(st_col, nrow(st_data)))
    is_data_row <- c(is_data_row, rep(TRUE, nrow(st_data)))
    
    # Empty row if not last
    if (i < length(levels_present)) {
      empty_row <- as.data.frame(matrix(NA, nrow=1, ncol=ncol(sub_df_base)))
      colnames(empty_row) <- colnames(sub_df_base)
      empty_row$Regulon <- ""
      final_df <- rbind(final_df, empty_row)
      row_colors <- c(row_colors, NA)
      is_data_row <- c(is_data_row, FALSE)
    }
  }
  
  # Headers
  writeData(wb, sheet_name, "Regulon", startCol=1, startRow=1)
  mergeCells(wb, sheet_name, cols = 1, rows = 1:2)
  for (st_label in names(state_header_pos)) {
    pos <- state_header_pos[[st_label]]
    writeData(wb, sheet_name, st_label, startCol=pos$middle, startRow=1)
    addStyle(wb, sheet_name, border_style, rows = 1:(nrow(final_df)+2), cols = pos$start, stack = TRUE)
  }
  writeData(wb, sheet_name, "scAtlas AUC", startCol = auc_sc_start, startRow = 1)
  writeData(wb, sheet_name, "PDO AUC",     startCol = auc_pdo_start, startRow = 1)
  addStyle(wb, sheet_name, sc_header_style, rows = 1, cols = auc_sc_start:auc_sc_end, gridExpand = TRUE)
  addStyle(wb, sheet_name, pdo_header_style, rows = 1, cols = auc_pdo_start:auc_pdo_end, gridExpand = TRUE)
  writeData(wb, sheet_name, t(display_colnames[-1]), startCol=2, startRow=2, colNames=FALSE)
  addStyle(wb, sheet_name, boldStyle, rows=1:2, cols=1:ncol(final_df), gridExpand=TRUE, stack=TRUE)
  
  # Write Data
  writeData(wb, sheet_name, final_df, startCol=1, startRow=3, colNames=FALSE)
  
  # Row-specific formatting for Column 1
  for (r in seq_len(nrow(final_df))) {
    if (is_data_row[r]) {
      st_font_style <- createStyle(textDecoration = "bold", fontName = "Consolas", fontColour = row_colors[r])
      addStyle(wb, sheet_name, st_font_style, rows = r + 2, cols = 1, stack = TRUE)
    }
  }
  
  addStyle(wb, sheet_name, numStyle, rows = 3:(nrow(final_df)+2), cols = 2:ncol(final_df), gridExpand = TRUE, stack=TRUE)
  addStyle(wb, sheet_name, sep_style, rows = 1:(nrow(final_df)+2), cols = sep_col_idx, stack = TRUE)
  
  conditionalFormatting(wb, sheet_name, cols = 2:(sep_col_idx-1), rows = 3:(nrow(final_df)+2), 
                        style = c("#1D4E89", "#F8F4EC", "#B22222"), rule = c(min_gap, 0, max_gap), type = "colourScale")
  conditionalFormatting(wb, sheet_name, cols = (sep_col_idx+1):ncol(final_df), rows = 3:(nrow(final_df)+2), 
                        style = c("#FFFFFF", "#FB8A8A", "#B22222"), rule = c(min_auc, (min_auc + max_auc)/2, max_auc), type = "colourScale")
  addStyle(wb, sheet_name, border_style, rows = 1:(nrow(final_df)+2), cols = auc_sc_start, stack = TRUE)
  addStyle(wb, sheet_name, border_style, rows = 1:(nrow(final_df)+2), cols = auc_pdo_start, stack = TRUE)
  addStyle(wb, sheet_name, border_style, rows = 1:(nrow(final_df)+2), cols = sep_col_idx, stack = TRUE)
  setColWidths(wb, sheet_name, cols = 1, widths = 25)
  setColWidths(wb, sheet_name, cols = 2:ncol(final_df), widths = 12)
  setColWidths(wb, sheet_name, cols = sep_col_idx, widths = 4)
  freezePane(wb, sheet_name, firstActiveRow = 3, firstActiveCol = 2)
}

build_top5_sheet(wb2, "scAtlas Top5", sc_hm_df)
build_top5_sheet(wb2, "PDO Top5", pdo_hm_df)
build_top5_sheet(wb2, "Combined Top5", comb_regs_df)

out_xlsx_top5 <- file.path(pdo_dir, "Auto_scRef_PDO_scenic_top5_markers.xlsx")
saveWorkbook(wb2, out_xlsx_top5, overwrite = TRUE)
message("Saved TOP 5 comparative Excel table to ", out_xlsx_top5)

