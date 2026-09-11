####################
# Analysis registry (authoritative override):
#   Status: active
#   Script: analysis/cell_states/Auto_PDO_scAtlas_scenic_comparison.R
#   Methodology: analysis/methodology/cell_states/Auto_PDO_scAtlas_scenic_comparison_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Preserves a superseded implementation or analysis tied to superseded
#     inputs. Do not use its outputs as current centred-MP/state inputs. The
#     original historical inputs, outputs, and method notes remain below.
####################

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

setwd("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs")

sc_dir <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/final_mp_scenic"
pdo_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/final_mp_scenic"

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
sc_defined_states <- c("Classic proliferation", "Squamous-to-intestinal", "Glandular-to-intestinal", "Stress-adaptive")
pdo_defined_states <- c("Classic proliferation", "Columnar-to-intestinal", "Glandular differentiation", "Stress-adaptive")

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

  if (type == "MP") {
    meta_comb <- bind_rows(meta1, meta2)
    valid_groups <- c("Classic proliferation", "Squamous-to-intestinal", "Columnar-to-intestinal", "Glandular-to-intestinal", "Glandular differentiation", "Stress-adaptive")
    meta_comb <- meta_comb %>% filter(mp_group %in% valid_groups)
    
    meta_comb <- meta_comb %>%
      mutate(ordered_group = case_when(
        mp_group == "Classic proliferation" ~ "Classic proliferation",
        mp_group %in% c("Squamous-to-intestinal", "Columnar-to-intestinal") ~ "Intestinal",
        mp_group %in% c("Glandular-to-intestinal", "Glandular differentiation") ~ "Glandular",
        mp_group == "Stress-adaptive" ~ "Stress-adaptive",
        TRUE ~ "Other"
      ))
      
    meta_comb$ordered_group <- factor(meta_comb$ordered_group, levels = c("Classic proliferation", "Intestinal", "Glandular", "Stress-adaptive"))
    meta_comb$Dataset <- factor(meta_comb$Dataset, levels = c("scAtlas", "PDO"))
    meta_comb <- meta_comb %>% arrange(ordered_group, Dataset, final_mp_label)
    
    mat1_sub <- mat1[com_regs, intersect(colnames(mat1), meta_comb$final_mp_label[meta_comb$Dataset=="scAtlas"]), drop=FALSE]
    mat2_sub <- mat2[com_regs, intersect(colnames(mat2), meta_comb$final_mp_label[meta_comb$Dataset=="PDO"]), drop=FALSE]
    
    mat_comb <- cbind(mat1_sub, mat2_sub)
    mat_comb <- mat_comb[, meta_comb$final_mp_label, drop=FALSE]
    
    group_cols <- c(
      "Classic proliferation" = "#E41A1C",
      "Squamous-to-intestinal" = "#4DAF4A",
      "Columnar-to-intestinal" = "#4DAF4A",
      "Glandular-to-intestinal" = "#FF7F00",
      "Glandular differentiation" = "#FF7F00",
      "Stress-adaptive" = "#984EA3"
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
    sc_states <- c("Classic proliferation", "Squamous-to-intestinal", "Glandular-to-intestinal", "Stress-adaptive")
    pdo_states <- c("Classic proliferation", "Columnar-to-intestinal", "Glandular differentiation", "Stress-adaptive")
    
    col_order <- c(
      sc_states[1], pdo_states[1],
      sc_states[2], pdo_states[2],
      sc_states[3], pdo_states[3],
      sc_states[4], pdo_states[4]
    )
    
    mat1_sub <- mat1[com_regs, sc_states, drop=FALSE]
    mat2_sub <- mat2[com_regs, pdo_states, drop=FALSE]
    
    mat_comb <- cbind(
      mat1_sub[, 1, drop=FALSE], mat2_sub[, 1, drop=FALSE],
      mat1_sub[, 2, drop=FALSE], mat2_sub[, 2, drop=FALSE],
      mat1_sub[, 3, drop=FALSE], mat2_sub[, 3, drop=FALSE],
      mat1_sub[, 4, drop=FALSE], mat2_sub[, 4, drop=FALSE]
    )
    
    meta_comb <- data.frame(
      State = colnames(mat_comb),
      Dataset = rep(c("scAtlas", "PDO"), 4)
    )

    state_cols <- c(
      "Classic proliferation" = "#E41A1C",
      "Squamous-to-intestinal" = "#4DAF4A",
      "Columnar-to-intestinal" = "#4DAF4A",
      "Glandular-to-intestinal" = "#FF7F00",
      "Glandular differentiation" = "#FF7F00",
      "Stress-adaptive" = "#984EA3"
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
  
  mat_scaled <- t(scale(t(mat_comb)))
  mat_scaled[!is.finite(mat_scaled)] <- 0
  
  col_fun <- colorRamp2(c(0, 1.25, 2.5), c("#FFFFFF", "#FB8A8A", "#B22222"))
  
  pdf(out_pdf, width = 18, height = 15, useDingbats = FALSE)
  draw(
    Heatmap(
      mat_scaled,
      name = "Scaled RSS",
      col = col_fun,
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
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

get_auc_matrix <- function(auc_obj, metadata, group_col) {
  auc_mat <- AUCell::getAUC(auc_obj)
  groups <- unique(metadata[[group_col]])
  groups <- groups[!is.na(groups) & groups != ""]
  res <- list()
  for (g in groups) {
    cells <- metadata$cell[metadata[[group_col]] == g]
    cells <- intersect(cells, colnames(auc_mat))
    if (length(cells) > 0) {
      res[[g]] <- rowMeans(as.matrix(auc_mat[, cells, drop=FALSE]))
    } else {
      res[[g]] <- rep(NA, nrow(auc_mat))
    }
  }
  res_mat <- do.call(cbind, res)
  rownames(res_mat) <- format_regulon_name(rownames(auc_mat))
  if (any(duplicated(rownames(res_mat)))) {
    rsum <- rowsum(res_mat, rownames(res_mat))
    rcount <- table(rownames(res_mat))
    res_mat <- rsum / as.numeric(rcount[rownames(rsum)])
  }
  return(res_mat)
}

run_auc_heatmap <- function(mat1, mat2, meta1, meta2, out_pdf, title, type="MP") {
  com_regs <- intersect(rownames(mat1), rownames(mat2))
  message(sprintf("Found %d common regulons for %s AUC", length(com_regs), type))

  if (type == "MP") {
    meta_comb <- bind_rows(meta1, meta2)
    valid_groups <- c("Classic proliferation", "Squamous-to-intestinal", "Columnar-to-intestinal", "Glandular-to-intestinal", "Glandular differentiation", "Stress-adaptive")
    meta_comb <- meta_comb %>% filter(mp_group %in% valid_groups)
    
    meta_comb <- meta_comb %>%
      mutate(ordered_group = case_when(
        mp_group == "Classic proliferation" ~ "Classic proliferation",
        mp_group %in% c("Squamous-to-intestinal", "Columnar-to-intestinal") ~ "Intestinal",
        mp_group %in% c("Glandular-to-intestinal", "Glandular differentiation") ~ "Glandular",
        mp_group == "Stress-adaptive" ~ "Stress-adaptive",
        TRUE ~ "Other"
      ))
      
    meta_comb$ordered_group <- factor(meta_comb$ordered_group, levels = c("Classic proliferation", "Intestinal", "Glandular", "Stress-adaptive"))
    meta_comb$Dataset <- factor(meta_comb$Dataset, levels = c("scAtlas", "PDO"))
    meta_comb <- meta_comb %>% arrange(ordered_group, Dataset, final_mp_label)
    
    mat1_sub <- mat1[com_regs, intersect(colnames(mat1), meta_comb$final_mp_label[meta_comb$Dataset=="scAtlas"]), drop=FALSE]
    mat2_sub <- mat2[com_regs, intersect(colnames(mat2), meta_comb$final_mp_label[meta_comb$Dataset=="PDO"]), drop=FALSE]
    
    mat_comb <- cbind(mat1_sub, mat2_sub)
    mat_comb <- mat_comb[, meta_comb$final_mp_label, drop=FALSE]
    
    group_cols <- c(
      "Classic proliferation" = "#E41A1C",
      "Squamous-to-intestinal" = "#4DAF4A",
      "Columnar-to-intestinal" = "#4DAF4A",
      "Glandular-to-intestinal" = "#FF7F00",
      "Glandular differentiation" = "#FF7F00",
      "Stress-adaptive" = "#984EA3"
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
    sc_states <- c("Classic proliferation", "Squamous-to-intestinal", "Glandular-to-intestinal", "Stress-adaptive")
    pdo_states <- c("Classic proliferation", "Columnar-to-intestinal", "Glandular differentiation", "Stress-adaptive")
    
    col_order <- c(
      sc_states[1], pdo_states[1],
      sc_states[2], pdo_states[2],
      sc_states[3], pdo_states[3],
      sc_states[4], pdo_states[4]
    )
    
    mat1_sub <- mat1[com_regs, sc_states, drop=FALSE]
    mat2_sub <- mat2[com_regs, pdo_states, drop=FALSE]
    
    mat_comb <- cbind(
      mat1_sub[, 1, drop=FALSE], mat2_sub[, 1, drop=FALSE],
      mat1_sub[, 2, drop=FALSE], mat2_sub[, 2, drop=FALSE],
      mat1_sub[, 3, drop=FALSE], mat2_sub[, 3, drop=FALSE],
      mat1_sub[, 4, drop=FALSE], mat2_sub[, 4, drop=FALSE]
    )
    
    meta_comb <- data.frame(
      State = colnames(mat_comb),
      Dataset = rep(c("scAtlas", "PDO"), 4)
    )

    state_cols <- c(
      "Classic proliferation" = "#E41A1C",
      "Squamous-to-intestinal" = "#4DAF4A",
      "Columnar-to-intestinal" = "#4DAF4A",
      "Glandular-to-intestinal" = "#FF7F00",
      "Glandular differentiation" = "#FF7F00",
      "Stress-adaptive" = "#984EA3"
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
  
  vals <- as.numeric(mat_comb)
  vals <- vals[!is.na(vals) & is.finite(vals)]
  q95 <- if (length(vals) > 0) quantile(vals, 0.95, na.rm=TRUE) else 0.1
  q95 <- max(q95, 0.06)
  
  col_fun_auc <- colorRamp2(c(0, 0.025, q95), c("#1D4E89", "#F8F4EC", "#B22222"))
  
  pdf(out_pdf, width = 18, height = 15, useDingbats = FALSE)
  draw(
    Heatmap(
      mat_comb,
      name = "AUCell Score",
      col = col_fun_auc,
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = FALSE,
      show_column_dend = FALSE,
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

message("Generating AUC matrices...")
sc_mp_auc <- get_auc_matrix(sc_auc_mat, sc_metadata, "final_mp_label")
pdo_mp_auc <- get_auc_matrix(pdo_auc_mat, pdo_metadata, "final_mp_label")

sc_st_auc <- get_auc_matrix(sc_auc_mat, sc_metadata, "final_state")
pdo_st_auc <- get_auc_matrix(pdo_auc_mat, pdo_metadata, "final_state")

# AUC Heatmap for MPs
run_auc_heatmap(
  sc_mp_auc, pdo_mp_auc, 
  sc_mp_anno, pdo_mp_anno, 
  file.path(out_dir_plot, "Auto_scenic_comparison_MP_AUC_heatmap.pdf"),
  "SCENIC Regulon AUCell Comparison (scAtlas vs PDO MPs)",
  type="MP"
)

# AUC Heatmap for States
run_auc_heatmap(
  sc_st_auc, pdo_st_auc, 
  NULL, NULL, 
  file.path(out_dir_plot, "Auto_scenic_comparison_State_AUC_heatmap.pdf"),
  "SCENIC Regulon AUCell Comparison (scAtlas vs PDO States)",
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
  "Classic proliferation" = list(sc = "Classic proliferation", pdo = "Classic proliferation"),
  "Intestinal" = list(sc = "Squamous-to-intestinal", pdo = "Columnar-to-intestinal"),
  "Glandular" = list(sc = "Glandular-to-intestinal", pdo = "Glandular differentiation"),
  "Stress-adaptive" = list(sc = "Stress-adaptive", pdo = "Stress-adaptive")
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
  "Classic proliferation" = list(sc = "Classic proliferation", pdo = "Classic proliferation"),
  "Intestinal" = list(sc = "Squamous-to-intestinal", pdo = "Columnar-to-intestinal"),
  "Glandular" = list(sc = "Glandular-to-intestinal", pdo = "Glandular differentiation"),
  "Stress-adaptive" = list(sc = "Stress-adaptive", pdo = "Stress-adaptive")
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



# --- Section B: Separator ---
overview_df$Sep <- ""
sep_col_idx <- current_col
current_col <- current_col + 1

# --- Section C: Grouped AUC (Marker style) ---
state_abbrev <- c(
  "Classic proliferation" = "ClassProlif",
  "Intestinal" = "Intestinal",
  "Glandular" = "Glandular",
  "Stress-adaptive" = "StressAdapt"
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
  "Classic proliferation",
  "Squamous-to-intestinal",
  "Glandular-to-intestinal",
  "Stress-adaptive"
)

pdo_state_order <- c(
  "Classic proliferation",
  "Columnar-to-intestinal",
  "Glandular differentiation",
  "Stress-adaptive"
)

sc_state_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Squamous-to-intestinal" = "#4DAF4A",
  "Glandular-to-intestinal" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "Cancer-cell immune mimicry" = "#377EB8"
)

pdo_state_cols <- c(
  "Classic proliferation" = "#E41A1C",
  "Columnar-to-intestinal" = "#4DAF4A",
  "Glandular differentiation" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "ECM-remodelling" = "#A65628",
  "Motile-cilia differentiation" = "#F781BF",
  "Hybrid" = "black",
  "Unresolved" = "grey80"
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



comb_regs_df <- do.call(rbind, comb_regs_list)
comb_regs_df <- comb_regs_df[order(comb_regs_df$Score, decreasing=TRUE), ]
comb_regs_df <- comb_regs_df[!duplicated(comb_regs_df$Regulon), ]
comb_regs_df$State <- factor(comb_regs_df$State, levels = sc_state_order)
comb_regs_df <- comb_regs_df[order(comb_regs_df$State), ]
comb_regs <- comb_regs_df$Regulon

mat_sc <- sc_st_gap[comb_regs, sc_state_order, drop=FALSE]
colnames(mat_sc) <- sc_state_order

mat_pdo <- matrix(NA, nrow=length(comb_regs), ncol=length(pdo_state_order))
rownames(mat_pdo) <- comb_regs
colnames(mat_pdo) <- pdo_state_order
regs_in_pdo <- comb_regs %in% rownames(pdo_st_gap)
mat_pdo[regs_in_pdo, ] <- pdo_st_gap[comb_regs[regs_in_pdo], pdo_state_order, drop=FALSE]

  colnames(mat_sc) <- paste0("sc_", sc_state_order)
  colnames(mat_pdo) <- paste0("pdo_", pdo_state_order)
  col_order <- c(
    sc_state_order[1], pdo_state_order[1],
    sc_state_order[2], pdo_state_order[2],
    sc_state_order[3], pdo_state_order[3],
    sc_state_order[4], pdo_state_order[4]
  )
  col_order_uniq <- c(
    paste0("sc_", sc_state_order[1]), paste0("pdo_", pdo_state_order[1]),
    paste0("sc_", sc_state_order[2]), paste0("pdo_", pdo_state_order[2]),
    paste0("sc_", sc_state_order[3]), paste0("pdo_", pdo_state_order[3]),
    paste0("sc_", sc_state_order[4]), paste0("pdo_", pdo_state_order[4])
  )
  
  comb_plot_raw <- cbind(mat_sc, mat_pdo)
  comb_plot <- comb_plot_raw[, col_order_uniq, drop=FALSE]
  colnames(comb_plot) <- paste0(rep(c("scAtlas::", "PDO::"), 4), col_order)

row_ann_comb <- rowAnnotation(
  State = comb_regs_df$State,
  col = list(State = sc_state_cols),
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm")
)

col_split_factors <- factor(
  rep(c("scAtlas", "PDO"), 4),
  levels = c("scAtlas", "PDO")
)

col_state_vec <- c(
  sc_state_order[1], pdo_state_order[1],
  sc_state_order[2], pdo_state_order[2],
  sc_state_order[3], pdo_state_order[3],
  sc_state_order[4], pdo_state_order[4]
)
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

# Heatmap 4: Combined Heatmap based on AUCell Gap
comb_auc_regs_list <- list()
sc_st_auc_gap <- calc_rss_gap(sc_st_auc, sc_state_order)
pdo_st_auc_gap <- calc_rss_gap(pdo_st_auc, pdo_state_order)

auc_thresh <- quantile(c(as.numeric(sc_st_auc), as.numeric(pdo_st_auc))[c(as.numeric(sc_st_auc), as.numeric(pdo_st_auc)) > 0], 0.05, na.rm=TRUE)

for (st_label in names(shared_states)) {
  sc_st <- shared_states[[st_label]]$sc
  pdo_st <- shared_states[[st_label]]$pdo
  
  common <- intersect(rownames(sc_st_auc_gap), rownames(pdo_st_auc_gap))
  sc_vals <- sc_st_auc_gap[common, sc_st]
  pdo_vals <- pdo_st_auc_gap[common, pdo_st]
  
  sc_auc_val <- sc_st_auc[common, sc_st]
  pdo_auc_val <- pdo_st_auc[common, pdo_st]
  
  has_support <- sc_auc_val > auc_thresh & pdo_auc_val > auc_thresh
  
  if (any(has_support)) {
    comb_score <- (sc_vals[has_support] + pdo_vals[has_support]) / 2
    top5 <- names(sort(comb_score, decreasing=TRUE)[1:min(5, length(comb_score))])
    comb_auc_regs_list[[sc_st]] <- data.frame(Regulon = top5, State = sc_st, Score = comb_score[top5], stringsAsFactors=FALSE)
  }
}

if (length(comb_auc_regs_list) > 0) {
  comb_auc_regs_df <- do.call(rbind, comb_auc_regs_list)
  comb_auc_regs_df <- comb_auc_regs_df[order(comb_auc_regs_df$Score, decreasing=TRUE), ]
  comb_auc_regs_df <- comb_auc_regs_df[!duplicated(comb_auc_regs_df$Regulon), ]
  comb_auc_regs_df$State <- factor(comb_auc_regs_df$State, levels = sc_state_order)
  comb_auc_regs_df <- comb_auc_regs_df[order(comb_auc_regs_df$State), ]
  comb_auc_regs <- comb_auc_regs_df$Regulon

  mat_sc_auc_sel <- sc_st_auc[comb_auc_regs, sc_state_order, drop=FALSE]
  colnames(mat_sc_auc_sel) <- sc_state_order

  mat_pdo_auc_sel <- matrix(NA, nrow=length(comb_auc_regs), ncol=length(pdo_state_order))
  rownames(mat_pdo_auc_sel) <- comb_auc_regs
  colnames(mat_pdo_auc_sel) <- pdo_state_order
  regs_in_pdo_auc <- comb_auc_regs %in% rownames(pdo_st_auc)
  mat_pdo_auc_sel[regs_in_pdo_auc, ] <- pdo_st_auc[comb_auc_regs[regs_in_pdo_auc], pdo_state_order, drop=FALSE]

  colnames(mat_sc_auc_sel) <- paste0("sc_", sc_state_order)
  colnames(mat_pdo_auc_sel) <- paste0("pdo_", pdo_state_order)
  comb_plot_raw_auc_sel <- cbind(mat_sc_auc_sel, mat_pdo_auc_sel)
  comb_plot_auc_sel <- comb_plot_raw_auc_sel[, col_order_uniq, drop=FALSE]
  colnames(comb_plot_auc_sel) <- paste0(rep(c("scAtlas::", "PDO::"), 4), col_order)

  row_ann_comb_auc <- rowAnnotation(
    State = comb_auc_regs_df$State,
    col = list(State = sc_state_cols),
    show_annotation_name = FALSE,
    simple_anno_size = unit(4, "mm")
  )

  vals_auc_p4 <- as.numeric(comb_plot_auc_sel)
  vals_auc_p4 <- vals_auc_p4[!is.na(vals_auc_p4) & is.finite(vals_auc_p4)]
  q95_p4 <- if (length(vals_auc_p4) > 0) quantile(vals_auc_p4, 0.95, na.rm=TRUE) else 0.1
  q95_p4 <- max(q95_p4, 0.06)
  
  col_fun_auc_p4 <- colorRamp2(c(0, 0.025, q95_p4), c("#1D4E89", "#F8F4EC", "#B22222"))

  ht_comb_auc <- Heatmap(
    comb_plot_auc_sel, name = "AUCell\nScore",
    top_annotation = top_ann_comb, left_annotation = row_ann_comb_auc, col = col_fun_auc_p4,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
    row_split = comb_auc_regs_df$State, row_title_rot = 0,
    column_split = col_split_factors,
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45, border = TRUE
  )
}

# Heatmap 5: Combined Heatmap based on Highest AUCell Activity (No Gap)
comb_auc_act_regs_list <- list()

for (st_label in names(shared_states)) {
  sc_st <- shared_states[[st_label]]$sc
  pdo_st <- shared_states[[st_label]]$pdo
  
  common <- intersect(rownames(sc_st_auc), rownames(pdo_st_auc))
  
  sc_auc_val <- sc_st_auc[common, sc_st]
  pdo_auc_val <- pdo_st_auc[common, pdo_st]
  
  comb_score <- (sc_auc_val + pdo_auc_val) / 2
  top5 <- names(sort(comb_score, decreasing=TRUE)[1:min(5, length(comb_score))])
  comb_auc_act_regs_list[[sc_st]] <- data.frame(Regulon = top5, State = sc_st, Score = comb_score[top5], stringsAsFactors=FALSE)
}

if (length(comb_auc_act_regs_list) > 0) {
  comb_auc_act_regs_df <- do.call(rbind, comb_auc_act_regs_list)
  comb_auc_act_regs_df <- comb_auc_act_regs_df[order(comb_auc_act_regs_df$Score, decreasing=TRUE), ]

  comb_auc_act_regs_df$State <- factor(comb_auc_act_regs_df$State, levels = sc_state_order)
  comb_auc_act_regs_df <- comb_auc_act_regs_df[order(comb_auc_act_regs_df$State), ]
  comb_auc_act_regs <- comb_auc_act_regs_df$Regulon

  mat_sc_auc_act_sel <- sc_st_auc[comb_auc_act_regs, sc_state_order, drop=FALSE]
  colnames(mat_sc_auc_act_sel) <- sc_state_order

  mat_pdo_auc_act_sel <- matrix(NA, nrow=length(comb_auc_act_regs), ncol=length(pdo_state_order))
  rownames(mat_pdo_auc_act_sel) <- comb_auc_act_regs
  colnames(mat_pdo_auc_act_sel) <- pdo_state_order
  regs_in_pdo_auc_act <- comb_auc_act_regs %in% rownames(pdo_st_auc)
  mat_pdo_auc_act_sel[regs_in_pdo_auc_act, ] <- pdo_st_auc[comb_auc_act_regs[regs_in_pdo_auc_act], pdo_state_order, drop=FALSE]

  colnames(mat_sc_auc_act_sel) <- paste0("sc_", sc_state_order)
  colnames(mat_pdo_auc_act_sel) <- paste0("pdo_", pdo_state_order)
  comb_plot_raw_auc_act_sel <- cbind(mat_sc_auc_act_sel, mat_pdo_auc_act_sel)
  comb_plot_auc_act_sel <- comb_plot_raw_auc_act_sel[, col_order_uniq, drop=FALSE]
  colnames(comb_plot_auc_act_sel) <- paste0(rep(c("scAtlas::", "PDO::"), 4), col_order)

  row_ann_comb_auc_act <- rowAnnotation(
    State = comb_auc_act_regs_df$State,
    col = list(State = sc_state_cols),
    show_annotation_name = FALSE,
    simple_anno_size = unit(4, "mm")
  )

  ht_comb_auc_act <- Heatmap(
    comb_plot_auc_act_sel, name = "AUCell\nScore\n(Activity Sel)",
    top_annotation = top_ann_comb, left_annotation = row_ann_comb_auc_act, col = col_fun_auc_p4,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
    row_split = comb_auc_act_regs_df$State, row_title_rot = 0,
    column_split = col_split_factors,
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 45, border = TRUE
  )
}

pdf(file.path(out_dir_plot, "Auto_scenic_comparison_RSS_heatmaps.pdf"), width = 17, height = 12, useDingbats = FALSE)
draw(ht_sc)
draw(ht_pdo)
draw(ht_comb)
if (length(comb_auc_regs_list) > 0) {
  draw(ht_comb_auc)
}
if (length(comb_auc_act_regs_list) > 0) {
  draw(ht_comb_auc_act)
}
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

####################
# Addendum: Independent state-vs-rest RSS concordance
####################
library(patchwork)
library(ggrepel)
library(tidyr)

message("Generating RSS Gap state-vs-state concordance correlation and scatter plots...")

com_regs_all <- intersect(rownames(sc_st_rss), rownames(pdo_st_rss))
sc_all_states <- colnames(sc_st_rss)
pdo_all_states <- colnames(pdo_st_rss)

sc_rss_mat <- sc_st_rss[com_regs_all, , drop=FALSE]
pdo_rss_mat <- pdo_st_rss[com_regs_all, , drop=FALSE]

calc_gap <- function(mat) {
  res <- mat
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      res[i,j] <- mat[i,j] - max(mat[i,-j], na.rm=TRUE)
    }
  }
  res
}

sc_gap_mat <- calc_gap(sc_rss_mat)
pdo_gap_mat <- calc_gap(pdo_rss_mat)

rss_cor <- cor(pdo_gap_mat, sc_gap_mat, method = "spearman", use = "pairwise.complete.obs")

expected_state_map <- c(
  "Classic proliferation" = "Classic proliferation",
  "Columnar-to-intestinal" = "Squamous-to-intestinal",
  "Glandular differentiation" = "Glandular-to-intestinal",
  "Stress-adaptive" = "Stress-adaptive"
)

# Heatmap
matrix_to_long <- function(mat, row_name, col_name, value_name) {
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c(row_name, col_name, value_name)
  df
}

plot_df <- matrix_to_long(rss_cor, "row_state", "col_state", "value")
plot_df$row_state <- factor(plot_df$row_state, levels = rev(pdo_all_states))
plot_df$col_state <- factor(plot_df$col_state, levels = sc_all_states)
plot_df$label <- sprintf("%.2f", plot_df$value)
plot_df$expected <- as.character(plot_df$row_state) %in% names(expected_state_map) &
  expected_state_map[as.character(plot_df$row_state)] == as.character(plot_df$col_state)

p_heatmap <- ggplot(plot_df, aes(x = col_state, y = row_state, fill = value)) +
  geom_tile(aes(color = expected), linewidth = 1.2) +
  geom_text(aes(label = label), size = 3.8, fontface = "bold") +
  scale_color_manual(values = c("FALSE" = "white", "TRUE" = "black"), guide = "none") +
  labs(
    title = "Independent state-vs-rest regulon specificity gap correlation",
    x = "scRef state", 
    y = "PDO state", 
    fill = "Spearman"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "right"
  ) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1))

# Scatters
make_rss_scatter <- function(pdo_state, sc_state) {
  df <- data.frame(
    Regulon = com_regs_all,
    pdo_gap = pdo_gap_mat[, pdo_state],
    sc_gap = sc_gap_mat[, sc_state],
    stringsAsFactors = FALSE
  )
  rho <- cor(df$pdo_gap, df$sc_gap, method = "spearman", use = "pairwise.complete.obs")
  
  df$rank_score <- df$pdo_gap + df$sc_gap
  df <- df[order(df$rank_score, decreasing = TRUE), ]
  
  df$label <- ""
  df$label[1:8] <- df$Regulon[1:8]
  
  max_val <- max(abs(c(df$pdo_gap, df$sc_gap)), na.rm=TRUE) * 1.05
  
  st_col <- pdo_state_cols[pdo_state]
  if(is.na(st_col)) st_col <- "grey50"
  
  ggplot(df, aes(x = sc_gap, y = pdo_gap)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point(color = st_col, size = 1, alpha = 0.4) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = Inf, min.segment.length = 0, seed = 1) +
    annotate("text", x = -0.96 * max_val, y = 0.96 * max_val, hjust = 0, vjust = 1, 
             label = sprintf("Spearman rho = %.2f", rho), size = 3, fontface = "bold") +
    coord_fixed(xlim = c(-max_val, max_val), ylim = c(-max_val, max_val), expand = FALSE) +
    labs(title = sc_state, x = "scRef RSS Specificity Gap", y = "PDO RSS Specificity Gap") +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
      plot.margin = margin(5, 4, 5, 4)
    )
}

pdf_out <- file.path(pdo_dir, "scAtlas_comparison", "Auto_scenic_comparison_concordance_scatter.pdf")
pdf(pdf_out, width = 16, height = 9, useDingbats = FALSE)
print(p_heatmap)

for (pdo_st in pdo_all_states) {
  plots <- lapply(sc_all_states, function(sc_st) make_rss_scatter(pdo_st, sc_st))
  if (length(plots) > 6) {
    chunk1 <- plots[1:5]
    chunk2 <- plots[6:length(plots)]
    print(wrap_plots(chunk1, nrow = 1) + plot_annotation(title = paste0("PDO state: ", pdo_st, " vs scAtlas (Part 1)"), theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))))
    print(wrap_plots(chunk2, nrow = 1) + plot_annotation(title = paste0("PDO state: ", pdo_st, " vs scAtlas (Part 2)"), theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))))
  } else {
    print(wrap_plots(plots, nrow = 1) + plot_annotation(title = paste0("PDO state: ", pdo_st, " vs all scAtlas states"), theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5))))
  }
}
dev.off()
message("Saved concordance plots to ", pdf_out)
