####################
# Analysis registry:
#   Status: active terminal survival analysis
#   Script: analysis/metaprograms/centred/Auto_08_tcga_mp_survival_volcano_centred.R
#   Methodology: analysis/methodology/metaprograms/centred/Auto_centred_ordering_and_state_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Scores current centred-refined PDO MP gene sets in TCGA ESCA bulk TPM and
#     estimates MP-specific overall-survival associations for EAC cases using
#     the Cox-model/split sensitivity analyses implemented below.
#
#   Inputs:
#     - scRef TCGA metadata: /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds
#     - scRef TCGA TPM: /rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt
#     - live: PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds
#
#   Outputs (live: PDOs_outs/centred_mp_refinement/survival/):
#     - Auto_task2_centred_mp_survival_volcano_centred_mps.pdf
#     - Auto_task2_centred_mp_survival_mp_cox_methods_splits.csv
#   Downstream use: none; terminal survival table and presentation figure.
#   Cache/replot behavior: rebuilds from persistent TCGA inputs and MP genes.
#   Run command: qsub Auto_centred_07_08.sh
#
#   Conda env: dmtcp
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(survival)
library(GSVA)

live_base <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs"
ephemeral_base <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs"
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")

task_prefix <- "task2_centred_mp"
out_dir <- file.path(live_base, "centred_mp_refinement", "survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

open_pdf_device <- function(path, width, height) {
  grDevices::cairo_pdf(filename = path, width = width, height = height, onefile = TRUE)
}

write_grob_pdf <- function(path, grob_list, width, height) {
  grob_list <- grob_list[!vapply(grob_list, is.null, logical(1))]
  if (length(grob_list) == 0) {
    message("No grobs available to write: ", path)
    return()
  }
  open_pdf_device(path, width = width, height = height)
  for (g in grob_list) {
    grid::grid.newpage()
    grid::grid.draw(g)
  }
  dev.off()
}

infer_histology <- function(type_vec) {
  t <- tolower(as.character(type_vec))
  out <- rep("Other", length(t))
  out[grepl("adeno", t)] <- "EAC"
  out
}

run_gsva <- function(expr_mat, gene_sets) {
  if (is.null(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) < 10) return(NULL)
  gs <- lapply(gene_sets, function(g) intersect(unique(g), rownames(expr_mat)))
  gs <- gs[sapply(gs, length) >= 5]
  if (length(gs) == 0) return(NULL)
  gsva(expr_mat, gs, method = "gsva", kcdf = "Gaussian")
}

run_cox <- function(df, feature_cols, mode_name, method_name, feature_type, split_method = "continuous") {
  out <- list()
  for (feat in feature_cols) {
    if (!feat %in% colnames(df)) next
    d <- df %>% filter(!is.na(OS_time), !is.na(OS_event), !is.na(.data[[feat]]), HistologyGroup == "EAC")
    if (nrow(d) < 20 || var(d[[feat]], na.rm = TRUE) == 0) next
    
    if (split_method == "median") {
      med_val <- median(d[[feat]], na.rm = TRUE)
      d$split_val <- factor(ifelse(d[[feat]] > med_val, "High", "Low"), levels = c("Low", "High"))
    } else if (split_method == "q1q4") {
      quants <- quantile(d[[feat]], probs = c(0.25, 0.75), na.rm = TRUE)
      d <- d %>% filter(.data[[feat]] <= quants[1] | .data[[feat]] >= quants[2])
      if (nrow(d) < 20) next
      d$split_val <- factor(ifelse(d[[feat]] >= quants[2], "High", "Low"), levels = c("Low", "High"))
    } else {
      d$split_val <- d[[feat]]
    }
    
    if (nrow(d) < 20 || var(as.numeric(d$split_val), na.rm = TRUE) == 0) next
    form <- if (split_method == "continuous") {
      as.formula(paste0("Surv(OS_time, OS_event) ~ `", feat, "`"))
    } else {
      as.formula("Surv(OS_time, OS_event) ~ split_val")
    }
    
    fit <- try(coxph(form, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) next
    ss <- summary(fit)
    out[[feat]] <- data.frame(
      mode = mode_name,
      method = method_name,
      cohort = "EAC",
      feature_type = feature_type,
      feature = feat,
      split_method = split_method,
      HR = ss$coefficients[1, "exp(coef)"],
      P_value = ss$coefficients[1, "Pr(>|z|)"],
      n = fit$n,
      events = fit$nevent,
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  bind_rows(out)
}

plot_volcano <- function(df, ttl) {
  if (nrow(df) == 0) return(NULL)
  pdat <- df %>%
    mutate(sig = P_value < 0.05, log2HR = log2(HR), neglog10 = -log10(P_value))
  ggplot(pdat, aes(log2HR, neglog10)) +
    geom_point(aes(color = sig), size = 2.8, alpha = 0.9) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey45") +
    geom_text_repel(aes(label = feature), size = 2.8, max.overlaps = 100) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(title = ttl, x = "log2(HR)", y = "-log10(p)")
}

make_tcga_page <- function(plot1, plot2, plot3, page_title) {
  gridExtra::arrangeGrob(
    plot1, plot2, plot3,
    ncol = 3,
    top = grid::textGrob(page_title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  )
}

make_feature_label <- function(x, feature_type, mp_desc) {
  if (feature_type == "MP") {
    d <- unname(mp_desc[x])
    d[is.na(d)] <- x[is.na(d)]
    return(paste0(x, " ", d))
  }
  x
}

mp_desc_map <- PDO_MP_DESCRIPTIONS

# Load Metadata
meta_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction/intermediate/Auto_tcga_esca_meta.rds"
if (!file.exists(meta_path)) {
  stop("Metadata file not found: ", meta_path)
}
meta_tcga <- readRDS(meta_path)
if (!"HistologyGroup" %in% colnames(meta_tcga)) {
  meta_tcga$HistologyGroup <- infer_histology(meta_tcga$type)
}

# Load Whole TCGA
whole_path <- "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/TCGA/esca_gdc_reconstruction/tables/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt"
if (!file.exists(whole_path)) {
  stop("Whole TPM file not found: ", whole_path)
}
tpm_df <- data.table::fread(whole_path)
tpm_whole <- as.matrix(tpm_df[, -1])
rownames(tpm_whole) <- tpm_df$GeneSymbol

mp_genes_path <- file.path(live_base, "centred_mp_refinement", "merged_refined_mp_genes.rds")
if (!file.exists(mp_genes_path)) {
  mp_genes_path <- file.path(ephemeral_base, "centred_mp_refinement", "intermediate", "merged_refined_mp_genes.rds")
}
if(!file.exists(mp_genes_path)) {
  stop("MP genes file not found: ", mp_genes_path)
}
mp.genes <- readRDS(mp_genes_path)
retained_mps <- names(mp.genes)

state_groups <- PDO_MP_STATE_GROUPS

state.genes <- lapply(state_groups, function(mps) {
  unique(unlist(mp.genes[intersect(mps, names(mp.genes))], use.names = FALSE))
})
state.genes <- state.genes[sapply(state.genes, length) >= 5]

all_res <- list()
split_methods <- c("continuous", "median", "q1q4")
panel_results_mp <- list()
panel_results_state <- list()

mode_name <- "noreg"
method_name <- "whole_tcga"
expr_mat <- tpm_whole

mp_gs <- run_gsva(expr_mat, mp.genes)
st_gs <- run_gsva(expr_mat, state.genes)

if (!is.null(mp_gs) || !is.null(st_gs)) {
  merged_df <- meta_tcga %>% filter(sample_type_code == "01")
  
  if (!is.null(mp_gs)) {
    mp_df <- as.data.frame(t(mp_gs))
    mp_df$sample_barcode <- rownames(mp_df)
    merged_df <- merged_df %>% left_join(mp_df, by = "sample_barcode")
    mp_cols <- intersect(colnames(as.data.frame(t(mp_gs))), colnames(merged_df))
  } else {
    mp_cols <- character(0)
  }
  
  if (!is.null(st_gs)) {
    st_df <- as.data.frame(t(st_gs))
    st_df$sample_barcode <- rownames(st_df)
    merged_df <- merged_df %>% left_join(st_df, by = "sample_barcode", suffix = c("", "_state"))
    st_cols <- intersect(colnames(as.data.frame(t(st_gs))), colnames(merged_df))
  } else {
    st_cols <- character(0)
  }
  
  for (sm in split_methods) {
    if (length(mp_cols) > 0) {
      cox_mp <- run_cox(merged_df, mp_cols, mode_name, method_name, "MP", split_method = sm)
      all_res[[paste("MP", sm, sep = "_")]] <- cox_mp
      
      this_mp <- cox_mp %>% filter(cohort == "EAC")
      if (nrow(this_mp) > 0) {
        mp_levels <- make_feature_label(retained_mps, "MP", mp_desc_map)
        this_mp$feature <- make_feature_label(this_mp$feature, "MP", mp_desc_map)
        all_levels <- unique(c(mp_levels, as.character(this_mp$feature)))
        this_mp <- this_mp %>% mutate(feature = factor(feature, levels = all_levels))
        panel_results_mp[[sm]] <- plot_volcano(
          this_mp,
          paste0("whole_tcga MP volcano (", sm, ")")
        )
      } else {
        panel_results_mp[[sm]] <- NULL
      }
    }
    
    if (length(st_cols) > 0) {
      cox_st <- run_cox(merged_df, st_cols, mode_name, method_name, "State", split_method = sm)
      all_res[[paste("State", sm, sep = "_")]] <- cox_st
      
      this_st <- cox_st %>% filter(cohort == "EAC")
      if (nrow(this_st) > 0) {
        this_st <- this_st %>% mutate(feature = factor(feature, levels = names(state_groups)))
        panel_results_state[[sm]] <- plot_volcano(
          this_st,
          paste0("whole_tcga State volcano (", sm, ")")
        )
      } else {
        panel_results_state[[sm]] <- NULL
      }
    }
  }
}

volcano_page_mp <- make_tcga_page(
  panel_results_mp[["continuous"]],
  panel_results_mp[["median"]],
  panel_results_mp[["q1q4"]],
  "Centred MP volcano: Whole TCGA"
)

volcano_page_state <- make_tcga_page(
  panel_results_state[["continuous"]],
  panel_results_state[["median"]],
  panel_results_state[["q1q4"]],
  "Centred State volcano: Whole TCGA"
)

write_grob_pdf(
  file.path(out_dir, paste0("Auto_", task_prefix, "_survival_volcano_centred_mps.pdf")),
  grob_list = list(volcano_page_mp, volcano_page_state),
  width = 18,
  height = 8
)

cox_res <- bind_rows(all_res)
if (nrow(cox_res) > 0) {
  cox_res$padj <- ave(
    cox_res$P_value,
    interaction(cox_res$mode, cox_res$method, cox_res$cohort, cox_res$feature_type, cox_res$split_method),
    FUN = function(x) p.adjust(x, method = "BH")
  )
  write.csv(
    cox_res,
    file.path(out_dir, paste0("Auto_", task_prefix, "_survival_mp_cox_methods_splits.csv")),
    row.names = FALSE
  )
}

cat("=== Auto_08 Complete ===\nSaved filtered TCGA whole-bulk MP volcano outputs for PDO centred MPs.\n")
