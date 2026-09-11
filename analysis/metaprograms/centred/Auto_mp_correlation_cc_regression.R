####################
# Analysis registry:
#   Status: active terminal figure
#   Script: analysis/metaprograms/centred/Auto_mp_correlation_cc_regression.R
#   Methodology: none (direct replotting with OLS regression variant)
#   Map: analysis/ANALYSIS_MAP.md
#
#   Description:
#     Produces a 2-page PDF comparing PDO MP correlations with and without
#     OLS regression of cell-cycle MP UCell scores from non-CC MPs.
#       Page 1: PDO intra-dataset Fisher-Z sample-averaged Spearman correlation.
#               Left = original scores, Right = CC-regressed non-CC + raw CC.
#       Page 2: scRef vs PDO cross-dataset Fisher-Z Spearman correlation.
#               Left = original scores, Right = both scRef & PDO CC-regressed.
#               Cell-cycle MPs are shown for both datasets in both panels.
#
#     Regression approach (from scRef legacy_state_definition_approach_b_reg_noreg.R):
#       X <- cbind(1, CC_scores);  B <- solve(X'X) %*% X'Y;  Y_reg <- Y - X %*% B
#
#   Inputs:
#     - live: PDOs_outs/centred_mp_refinement/merged_refined_ucell_scores.rds
#     - live: PDOs_outs/PDOs_merged.rds
#     - live: PDOs_outs/Auto_PDO_ucell_scATLAS_MPs_in_PDOs.rds
#     - live: scRef_Pipeline/analysis/shared/scRef_config.R
#     - live: PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R
#
#   Outputs (live: PDOs_outs/centred_mp_refinement/figures/):
#     Auto_mp_correlation_cc_regression.pdf   (2 pages)
#     Auto_mp_correlation_cc_regression.png   (page 1 only)
#
#   Downstream use: terminal figure; no downstream scripts depend on this.
#   Cache/replot behavior: reads existing UCell caches; never reruns NMF or scoring.
#   Run command: qsub Auto_replot_cc_regression.sh
#   Conda env: dmtcp
####################

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(grid)
})

source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")
source("/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/shared/scRef_config.R")

# === Paths ===
live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
fig_dir <- file.path(live_base, "centred_mp_refinement", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ucell_file <- file.path(live_base, "centred_mp_refinement", "merged_refined_ucell_scores.rds")
seurat_file <- file.path(live_base, "PDOs_merged.rds")
sc_ucell_in_pdo_file <- file.path(live_base, "Auto_PDO_ucell_scATLAS_MPs_in_PDOs.rds")

required_inputs <- c(ucell_file, seurat_file, sc_ucell_in_pdo_file)
missing <- required_inputs[!file.exists(required_inputs)]
if (length(missing) > 0) stop("Missing required input(s): ", paste(missing, collapse = ", "))

# === Load data ===
cat("Loading inputs...\n")
merged_ucell <- readRDS(ucell_file)
pdos_merged <- readRDS(seurat_file)
sc_ucell_raw <- readRDS(sc_ucell_in_pdo_file)  # MP x cells

# === PDO MP definitions from shared config ===
cc_mps <- PDO_CELL_CYCLE_MPS  # c("MP11", "MP1", "MP2", "MP3")
state_groups_with_cc <- c(list("Cell cycle" = cc_mps), PDO_MP_STATE_GROUPS)
strict_order <- unlist(state_groups_with_cc, use.names = FALSE)
strict_order <- strict_order[strict_order %in% colnames(merged_ucell)]
non_cc_mps <- setdiff(strict_order, cc_mps)

mp_to_state <- unlist(lapply(names(state_groups_with_cc), function(s) {
  setNames(rep(s, length(state_groups_with_cc[[s]])), state_groups_with_cc[[s]])
}))

state_cols_pdo <- c("Cell cycle" = "#6B7280", PDO_STATE_COLORS[PDO_STATE_ORDER])

mp_desc_map <- setNames(
  paste(names(PDO_MP_DESCRIPTIONS), PDO_MP_DESCRIPTIONS, sep = ": "),
  names(PDO_MP_DESCRIPTIONS)
)

display_label <- function(x) {
  desc <- unname(mp_desc_map[x])
  desc[is.na(desc)] <- x[is.na(desc)]
  ifelse(desc == x, x, desc)
}

# === scRef MP definitions from shared config ===
sc_cc_mps <- SCREF_CC_MPS  # c("MP1", "MP5", "MP13+")
sc_state_groups_with_cc <- c(list("Cell cycle" = sc_cc_mps), SCREF_STATE_GROUPS)
sc_strict_order <- unlist(sc_state_groups_with_cc, use.names = FALSE)

sc_mp_to_state <- unlist(lapply(names(sc_state_groups_with_cc), function(s) {
  setNames(rep(s, length(sc_state_groups_with_cc[[s]])), sc_state_groups_with_cc[[s]])
}))

sc_state_cols <- c("Cell cycle" = "#6B7280", SCREF_STATE_COLOURS[SCREF_PRIMARY_STATE_ORDER])

sc_mp_desc_map <- setNames(
  paste(names(SCREF_MP_DESCRIPTIONS), SCREF_MP_DESCRIPTIONS, sep = ": "),
  names(SCREF_MP_DESCRIPTIONS)
)

sc_display_label <- function(x) {
  desc <- unname(sc_mp_desc_map[x])
  desc[is.na(desc)] <- x[is.na(desc)]
  ifelse(desc == x, x, desc)
}

# === Sample vector ===
cell_meta <- pdos_merged@meta.data
cell_ids <- rownames(merged_ucell)
sample_vec <- cell_meta$orig.ident[match(cell_ids, rownames(cell_meta))]
sample_vec[is.na(sample_vec)] <- cell_ids[is.na(sample_vec)]

# === OLS regression function ===
regress_cc <- function(ucell_mat, cc_cols, non_cc_cols) {
  # Regress CC MP scores out of non-CC MPs using OLS.
  # Returns a matrix with same dims as ucell_mat but non_cc_cols replaced by residuals
  # and cc_cols unchanged.
  X_cc <- ucell_mat[, cc_cols, drop = FALSE]
  Y_other <- ucell_mat[, non_cc_cols, drop = FALSE]

  X <- cbind(Intercept = 1, X_cc)
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) {
    # Fallback: pseudoinverse via SVD
    svd_x <- svd(crossprod(X))
    svd_x$v %*% diag(1 / pmax(svd_x$d, 1e-10)) %*% t(svd_x$u)
  })
  B <- XtX_inv %*% crossprod(X, Y_other)
  Y_hat <- X %*% B
  Y_reg <- Y_other - Y_hat

  out <- ucell_mat
  out[, non_cc_cols] <- Y_reg
  out
}

# === Fisher-Z sample-averaged Spearman correlation ===
compute_fisher_cor <- function(score_mat, sample_vec, min_cells = 10) {
  score_mat <- as.matrix(score_mat)
  feature_names <- colnames(score_mat)
  samples <- unique(sample_vec)
  n_f <- length(feature_names)

  cor_array <- array(
    NA_real_,
    dim = c(n_f, n_f, length(samples)),
    dimnames = list(feature_names, feature_names, samples)
  )

  for (samp in samples) {
    idx <- which(sample_vec == samp)
    if (length(idx) < min_cells) next
    cor_array[, , samp] <- cor(score_mat[idx, , drop = FALSE], method = "spearman")
  }

  z_array <- atanh(pmin(pmax(cor_array, -0.999), 0.999))
  mean_rho <- matrix(NA_real_, n_f, n_f, dimnames = list(feature_names, feature_names))
  p_vals <- matrix(NA_real_, n_f, n_f, dimnames = list(feature_names, feature_names))

  for (i in seq_len(n_f)) {
    for (j in seq_len(n_f)) {
      if (i == j) {
        mean_rho[i, j] <- 1
        p_vals[i, j] <- 0
        next
      }
      zs <- z_array[i, j, ]
      zs <- zs[is.finite(zs)]
      if (length(zs) < 3) next
      mean_rho[i, j] <- tanh(mean(zs))
      tt <- tryCatch(t.test(zs), error = function(e) NULL)
      p_vals[i, j] <- if (!is.null(tt)) tt$p.value else NA_real_
    }
  }

  n_samples_used <- sum(vapply(samples, function(s) {
    length(which(sample_vec == s)) >= min_cells
  }, logical(1)))

  list(mean_rho = mean_rho, p_values = p_vals, n_samples = n_samples_used)
}

# === Cross-dataset Fisher-Z correlation (scRef rows x PDO cols) ===
compute_cross_fisher_cor <- function(pdo_mat, sc_mat, sample_vec, min_cells = 10) {
  # pdo_mat: cells x PDO MPs
  # sc_mat:  cells x scRef MPs (scored in PDO cells)
  # Returns mean_rho and p_vals: scRef (rows) x PDO (cols)
  common_cells <- intersect(rownames(pdo_mat), rownames(sc_mat))
  if (length(common_cells) < min_cells) stop("Fewer than min_cells common cells for cross-correlation.")

  pdo_mat <- as.matrix(pdo_mat[common_cells, , drop = FALSE])
  sc_mat <- as.matrix(sc_mat[common_cells, , drop = FALSE])
  sv <- sample_vec[common_cells]

  samples <- unique(sv)
  n_sc <- ncol(sc_mat)
  n_pdo <- ncol(pdo_mat)

  cor_list <- list()
  for (smp in samples) {
    idx <- which(sv == smp)
    if (length(idx) < min_cells) next
    cors <- matrix(NA_real_, n_sc, n_pdo)
    for (i in seq_len(n_sc)) {
      for (j in seq_len(n_pdo)) {
        if (sd(sc_mat[idx, i]) > 0 && sd(pdo_mat[idx, j]) > 0) {
          cors[i, j] <- cor(sc_mat[idx, i], pdo_mat[idx, j], method = "spearman")
        }
      }
    }
    rownames(cors) <- colnames(sc_mat)
    colnames(cors) <- colnames(pdo_mat)
    cor_list[[smp]] <- cors
  }

  if (length(cor_list) == 0) stop("No samples with >= min_cells for cross-correlation.")

  # Fisher-Z average
  sum_z <- matrix(0, n_sc, n_pdo)
  count_z <- matrix(0, n_sc, n_pdo)
  z_list <- list()
  for (smp in names(cor_list)) {
    z <- atanh(pmin(pmax(cor_list[[smp]], -0.999), 0.999))
    z_list[[smp]] <- z
    valid <- !is.na(z)
    sum_z[valid] <- sum_z[valid] + z[valid]
    count_z[valid] <- count_z[valid] + 1
  }
  mean_z <- sum_z / count_z
  mean_rho <- tanh(mean_z)

  p_vals <- matrix(1, n_sc, n_pdo)
  for (i in seq_len(n_sc)) {
    for (j in seq_len(n_pdo)) {
      vals <- vapply(z_list, function(z) z[i, j], numeric(1))
      vals <- vals[!is.na(vals)]
      if (length(vals) >= 3) {
        if (sd(vals) == 0) {
          p_vals[i, j] <- if (mean(vals) == 0) 1 else 1e-16
        } else {
          p_vals[i, j] <- t.test(vals, mu = 0)$p.value
        }
      }
    }
  }

  dimnames(mean_rho) <- list(colnames(sc_mat), colnames(pdo_mat))
  dimnames(p_vals) <- dimnames(mean_rho)

  list(mean_rho = mean_rho, p_values = p_vals, n_samples = length(cor_list))
}

# === Grab heatmap as grob (avoids viewport mixing with ComplexHeatmap) ===
grab_cor_heatmap <- function(mean_rho, p_vals, row_state_vec, col_state_vec,
                             row_state_cols, col_state_cols,
                             row_labels, col_labels,
                             title_text, n_samples,
                             hm_width = unit(9, "inch"),
                             hm_height = unit(9, "inch"),
                             row_legend_title = "State",
                             col_legend_title = "State",
                             show_row_legend = TRUE,
                             show_col_legend = TRUE) {

  plot_rho <- mean_rho
  rownames(plot_rho) <- row_labels[rownames(plot_rho)]
  colnames(plot_rho) <- col_labels[colnames(plot_rho)]

  col_cor <- colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

  ha_left <- rowAnnotation(
    State = row_state_vec,
    col = list(State = row_state_cols),
    show_annotation_name = FALSE,
    show_legend = show_row_legend,
    annotation_legend_param = list(title = row_legend_title)
  )
  ha_top <- HeatmapAnnotation(
    State = col_state_vec,
    col = list(State = col_state_cols),
    show_annotation_name = FALSE,
    show_legend = show_col_legend,
    annotation_legend_param = list(title = col_legend_title)
  )

  local_mean_rho <- mean_rho
  local_p_vals <- p_vals

  ht <- Heatmap(
    plot_rho,
    name = paste0("Mean Rho\n(", n_samples, " Samples)"),
    col = col_cor,
    rect_gp = gpar(col = "white", lwd = 1),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = ha_left,
    top_annotation = ha_top,
    row_split = row_state_vec,
    column_split = col_state_vec,
    column_title_rot = 20,
    column_title_side = "top",
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_title = NULL,
    row_names_side = "right",
    column_names_side = "bottom",
    column_names_rot = 30,
    row_names_gp = gpar(fontsize = 8, fontface = "bold"),
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),
    row_names_max_width = unit(100, "mm"),
    column_names_max_height = unit(100, "mm"),
    width = hm_width,
    height = hm_height,
    cell_fun = function(j, i, x, y, width, height, fill) {
      p <- local_p_vals[i, j]
      rho <- local_mean_rho[i, j]
      if (is.na(p) || is.na(rho)) {
        grid.text("NA", x, y, gp = gpar(fontsize = 6.5, col = "grey50"))
      } else if (p < 0.001) {
        grid.text(paste0(round(rho, 2), "\n***"), x, y, gp = gpar(fontsize = 6.5))
      } else if (p < 0.01) {
        grid.text(paste0(round(rho, 2), "\n**"), x, y, gp = gpar(fontsize = 6.5))
      } else if (p < 0.05) {
        grid.text(paste0(round(rho, 2), "\n*"), x, y, gp = gpar(fontsize = 6.5))
      } else {
        grid.text(round(rho, 2), x, y, gp = gpar(fontsize = 6.5))
      }
    },
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10)
    )
  )

  # Capture into a self-contained grob
  grid.grabExpr(
    draw(
      ht,
      heatmap_legend_side = "left",
      padding = unit(c(15, 15, 15, 15), "mm"),
      column_title = title_text,
      column_title_gp = gpar(fontsize = 14, fontface = "bold")
    ),
    wrap = TRUE
  )
}

# =====================================================================
# PAGE 1: PDO intra-dataset MP correlation - noreg vs reg
# =====================================================================
cat("Computing Page 1: PDO intra-dataset correlation...\n")

ucell_mat <- as.matrix(merged_ucell[, strict_order, drop = FALSE])

# --- No regression ---
cor_noreg <- compute_fisher_cor(ucell_mat, sample_vec)

# --- With regression ---
cc_in_data <- intersect(cc_mps, colnames(ucell_mat))
non_cc_in_data <- intersect(non_cc_mps, colnames(ucell_mat))
ucell_reg <- regress_cc(ucell_mat, cc_in_data, non_cc_in_data)
cor_reg <- compute_fisher_cor(ucell_reg, sample_vec)

# State vectors for annotation
state_vec_all <- factor(mp_to_state[strict_order], levels = names(state_groups_with_cc))
names_display_all <- setNames(display_label(strict_order), strict_order)

cat("Grabbing Page 1 heatmap grobs...\n")
grob_p1_left <- grab_cor_heatmap(
  cor_noreg$mean_rho, cor_noreg$p_values,
  row_state_vec = state_vec_all,
  col_state_vec = state_vec_all,
  row_state_cols = state_cols_pdo,
  col_state_cols = state_cols_pdo,
  row_labels = names_display_all,
  col_labels = names_display_all,
  title_text = "PDO MP Correlation - No Regression",
  n_samples = cor_noreg$n_samples,
  show_row_legend = TRUE,
  show_col_legend = FALSE
)

grob_p1_right <- grab_cor_heatmap(
  cor_reg$mean_rho, cor_reg$p_values,
  row_state_vec = state_vec_all,
  col_state_vec = state_vec_all,
  row_state_cols = state_cols_pdo,
  col_state_cols = state_cols_pdo,
  row_labels = names_display_all,
  col_labels = names_display_all,
  title_text = "PDO MP Correlation - CC Regressed",
  n_samples = cor_reg$n_samples,
  show_row_legend = FALSE,
  show_col_legend = FALSE
)

# =====================================================================
# PAGE 2: scRef vs PDO cross-dataset - noreg vs reg
# =====================================================================
cat("Computing Page 2: scRef vs PDO cross-dataset correlation...\n")

# Prepare scRef UCell scores in PDO cells (transpose: MP x cells -> cells x MP)
sc_mat_in_pdo <- t(as.matrix(sc_ucell_raw))

# Map scATLAS_ prefixed names to bare MP names
sc_bare_names <- sub("^scATLAS_", "", colnames(sc_mat_in_pdo))
colnames(sc_mat_in_pdo) <- sc_bare_names

# Keep only scRef MPs that are in the strict order
sc_keep <- intersect(sc_strict_order, colnames(sc_mat_in_pdo))
sc_mat_in_pdo <- sc_mat_in_pdo[, sc_keep, drop = FALSE]

# Prepare PDO UCell (cells x MP)
pdo_mat_for_cross <- as.matrix(merged_ucell[, strict_order, drop = FALSE])

# Common cells
common_cells <- intersect(rownames(pdo_mat_for_cross), rownames(sc_mat_in_pdo))
cat("  Cross-data common cells:", length(common_cells), "\n")

pdo_mat_for_cross <- pdo_mat_for_cross[common_cells, , drop = FALSE]
sc_mat_in_pdo <- sc_mat_in_pdo[common_cells, , drop = FALSE]

# Sample vector for common cells
sv_cross <- sample_vec[match(common_cells, cell_ids)]
names(sv_cross) <- common_cells

# --- No regression (left panel) ---
cross_noreg <- compute_cross_fisher_cor(pdo_mat_for_cross, sc_mat_in_pdo, sv_cross)

# --- With regression on both PDO and scRef sides (right panel) ---
pdo_mat_reg <- regress_cc(pdo_mat_for_cross, cc_in_data, non_cc_in_data)
sc_cc_in_data <- intersect(sc_cc_mps, colnames(sc_mat_in_pdo))
sc_non_cc_in_data <- setdiff(colnames(sc_mat_in_pdo), sc_cc_in_data)
sc_mat_reg <- regress_cc(sc_mat_in_pdo, sc_cc_in_data, sc_non_cc_in_data)

cross_reg <- compute_cross_fisher_cor(pdo_mat_reg, sc_mat_reg, sv_cross)

# State annotation vectors
state_vec_pdo_all <- factor(mp_to_state[strict_order], levels = names(state_groups_with_cc))
state_vec_sc <- factor(sc_mp_to_state[sc_keep], levels = names(sc_state_groups_with_cc))

# Display labels
pdo_labels_all <- setNames(display_label(strict_order), strict_order)
sc_labels <- setNames(sc_display_label(sc_keep), sc_keep)

cat("Grabbing Page 2 heatmap grobs...\n")
grob_p2_left <- grab_cor_heatmap(
  cross_noreg$mean_rho, cross_noreg$p_values,
  row_state_vec = state_vec_sc,
  col_state_vec = state_vec_pdo_all,
  row_state_cols = sc_state_cols,
  col_state_cols = state_cols_pdo,
  row_labels = sc_labels,
  col_labels = pdo_labels_all,
  title_text = "scRef vs PDO - No Regression",
  n_samples = cross_noreg$n_samples,
  row_legend_title = "scATLAS State",
  col_legend_title = "PDO State",
  show_row_legend = TRUE,
  show_col_legend = TRUE
)

grob_p2_right <- grab_cor_heatmap(
  cross_reg$mean_rho, cross_reg$p_values,
  row_state_vec = state_vec_sc,
  col_state_vec = state_vec_pdo_all,
  row_state_cols = sc_state_cols,
  col_state_cols = state_cols_pdo,
  row_labels = sc_labels,
  col_labels = pdo_labels_all,
  title_text = "scRef vs PDO - Both CC Regressed",
  n_samples = cross_reg$n_samples,
  row_legend_title = "scATLAS State",
  col_legend_title = "PDO State",
  show_row_legend = FALSE,
  show_col_legend = FALSE
)

# =====================================================================
# RENDER PDF & PNGs
# =====================================================================
pdf_path <- file.path(fig_dir, "Auto_mp_correlation_cc_regression.pdf")
png_p1_path <- file.path(fig_dir, "Auto_mp_correlation_cc_regression_page1.png")
png_p2_path <- file.path(fig_dir, "Auto_mp_correlation_cc_regression_page2.png")
png_default_path <- file.path(fig_dir, "Auto_mp_correlation_cc_regression.png")

cat("Drawing 2-page PDF...\n")

pdf(pdf_path, width = 38, height = 18, useDingbats = FALSE, onefile = TRUE)

# Page 1: PDO intra-dataset
grid.newpage()
pushViewport(viewport(x = 0, width = 0.5, just = "left"))
grid.draw(grob_p1_left)
popViewport()
pushViewport(viewport(x = 0.5, width = 0.5, just = "left"))
grid.draw(grob_p1_right)
popViewport()

# Page 2: scRef vs PDO cross-dataset
grid.newpage()
pushViewport(viewport(x = 0, width = 0.5, just = "left"))
grid.draw(grob_p2_left)
popViewport()
pushViewport(viewport(x = 0.5, width = 0.5, just = "left"))
grid.draw(grob_p2_right)
popViewport()

dev.off()

cat("Drawing Page 1 PNG...\n")
png(png_default_path, width = 38, height = 18, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(x = 0, width = 0.5, just = "left"))
grid.draw(grob_p1_left)
popViewport()
pushViewport(viewport(x = 0.5, width = 0.5, just = "left"))
grid.draw(grob_p1_right)
popViewport()
dev.off()

file.copy(png_default_path, png_p1_path, overwrite = TRUE)

cat("Drawing Page 2 PNG...\n")
png(png_p2_path, width = 38, height = 18, units = "in", res = 300, bg = "white")
grid.newpage()
pushViewport(viewport(x = 0, width = 0.5, just = "left"))
grid.draw(grob_p2_left)
popViewport()
pushViewport(viewport(x = 0.5, width = 0.5, just = "left"))
grid.draw(grob_p2_right)
popViewport()
dev.off()

cat("Saved:", pdf_path, "\n")
cat("Saved:", png_default_path, "\n")
cat("Saved:", png_p1_path, "\n")
cat("Saved:", png_p2_path, "\n")
cat("\n=== Auto CC Regression Correlation Heatmaps Complete ===\n")
