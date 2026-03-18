####################
# Adapted from scRef_Pipeline for PDOs
# Auto_survival_clinical_mps.R
####################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(survival)
library(survminer)
library(GSVA)

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

cat("=== Auto survival clinical MPs (PDO adapted) ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")

####################
# 1) Load MP definitions and UCell scores
####################
geneNMF.metaprograms <- readRDS("MP_outs_default.rds") # UPDATE: Change to optimal nMP result
ucell_scores <- readRDS("UCell_scores_filtered.rds") # UPDATE: Verify filename
tmdata_all <- readRDS("PDOs_merged.rds")

mp.genes <- geneNMF.metaprograms$metaprograms.genes
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
if (length(bad_mps) > 0) {
  mp.genes <- mp.genes[!names(mp.genes) %in% paste0("MP", bad_mps)]
}
retained_mps <- names(mp.genes)
cat("Retained MPs:", paste(retained_mps, collapse = ", "), "\n")

####################
# UPDATE AFTER geneNMF: Fill in MP descriptions and state groups
####################
mp_descriptions <- c(
  "MP6"  = "G2M Cell Cycle",
  "MP7"  = "DNA repair",
  "MP5"  = "MYC-related Proliferation",
  "MP1"  = "G2M checkpoint",
  "MP3"  = "G1S Cell Cycle",
  "MP8"  = "Columnar Progenitor",
  "MP10" = "Inflammatory Stress Epi.",
  "MP9"  = "ECM Remodeling Epi.",
  "MP4"  = "Intestinal Metaplasia"
)

state_groups <- list(
  "Classic Proliferative" = c("MP5"),
  "SMG-like Metaplasia"   = c("MP8"),
  "Stress-adaptive"       = c("MP10", "MP9"),
  "Basal to Intest. Meta" = c("MP4")
)

cc_mps <- c("MP6", "MP7", "MP1", "MP3")

tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
valid_cluster_ids <- as.numeric(gsub("\\D", "", retained_mps))
mp_tree_order <- unique(ordered_clusters)
mp_tree_order <- mp_tree_order[!is.na(mp_tree_order) & mp_tree_order %in% valid_cluster_ids]
mp_tree_order_names <- paste0("MP", mp_tree_order)
mp_cols <- intersect(mp_tree_order_names, colnames(ucell_scores))

####################
# 2) PDO clinical placeholder section (replacing TCGA loading)
####################
clinical_sheet <- readxl::read_excel(
  "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx"
)
clinical_df <- as.data.frame(clinical_sheet)
if (!"SUR" %in% colnames(clinical_df)) {
  rn <- rownames(clinical_df)
  if (!is.null(rn) && all(nchar(rn) > 0)) {
    clinical_df$SUR <- ifelse(grepl("^SUR", rn), rn, paste0("SUR", rn))
  }
}
clinical_df$SUR <- as.character(clinical_df$SUR)
clinical_df <- clinical_df %>% rename_with(~ paste0("clinical_", .x), -SUR)

meta_df <- tmdata_all@meta.data
meta_df$cell <- rownames(meta_df)
meta_df <- meta_df %>%
  mutate(
    orig.ident = as.character(orig.ident),
    SUR = ifelse("SUR" %in% colnames(.), as.character(SUR), sub("^(SUR[0-9]+).*", "\\1", as.character(orig.ident))),
    Batch = as.character(Batch)
  ) %>%
  left_join(clinical_df, by = "SUR") %>%
  filter(orig.ident != "SUR843T3_PDO")

response_vec <- if ("Clinical.response.at.OG.MDT..responder.non.responder" %in% colnames(meta_df)) {
  as.character(meta_df$Clinical.response.at.OG.MDT..responder.non.responder)
} else if ("clinical_Clinical.response.at.OG.MDT..responder.non.responder" %in% colnames(meta_df)) {
  as.character(meta_df$clinical_Clinical.response.at.OG.MDT..responder.non.responder)
} else {
  rep(NA_character_, nrow(meta_df))
}

meta_df$Clinical.response.at.OG.MDT..responder.non.responder <- response_vec

cell_ids <- Reduce(intersect, list(rownames(ucell_scores), rownames(meta_df)))
cell_long <- as.data.frame(ucell_scores[cell_ids, mp_cols, drop = FALSE]) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(mp_cols), names_to = "MP", values_to = "score") %>%
  left_join(meta_df %>% dplyr::select(cell, orig.ident, Batch, Clinical.response.at.OG.MDT..responder.non.responder), by = "cell") %>%
  mutate(
    Clinical.response.at.OG.MDT..responder.non.responder = ifelse(
      Clinical.response.at.OG.MDT..responder.non.responder %in% c("R", "Responder", "responder"),
      "Responder",
      ifelse(
        Clinical.response.at.OG.MDT..responder.non.responder %in% c("NR", "Nonresponder", "non-responder", "nonresponder"),
        "Nonresponder",
        Clinical.response.at.OG.MDT..responder.non.responder
      )
    )
  )

####################
# 3) Response-associated MP analysis (replacing Cox/KM)
####################
sample_mp <- cell_long %>%
  filter(!is.na(Clinical.response.at.OG.MDT..responder.non.responder), !is.na(score)) %>%
  group_by(orig.ident, Batch, Clinical.response.at.OG.MDT..responder.non.responder, MP) %>%
  summarise(mp_mean = mean(score, na.rm = TRUE), .groups = "drop")

run_wilcoxon <- function(df, mp_col = "MP", score_col = "mp_mean", group_col = "Clinical.response.at.OG.MDT..responder.non.responder") {
  out <- list()
  for (mp in unique(df[[mp_col]])) {
    d <- df %>% filter(.data[[mp_col]] == mp, !is.na(.data[[group_col]]), !is.na(.data[[score_col]]))
    if (nrow(d) < 4) next
    if (length(unique(d[[group_col]])) != 2) next
    gtab <- table(d[[group_col]])
    if (any(gtab < 2)) next

    wt <- suppressWarnings(wilcox.test(d[[score_col]] ~ d[[group_col]], exact = FALSE))
    grp_levels <- sort(unique(as.character(d[[group_col]])))
    ref_grp <- grp_levels[1]
    alt_grp <- grp_levels[2]
    eff <- median(d[[score_col]][d[[group_col]] == alt_grp], na.rm = TRUE) -
      median(d[[score_col]][d[[group_col]] == ref_grp], na.rm = TRUE)

    out[[mp]] <- data.frame(
      feature = mp,
      p_value = wt$p.value,
      effect = eff,
      ref_group = ref_grp,
      alt_group = alt_grp,
      n = nrow(d),
      stringsAsFactors = FALSE
    )
  }
  if (length(out) == 0) return(data.frame())
  res <- bind_rows(out)
  res$padj <- p.adjust(res$p_value, method = "BH")
  res
}

wilcox_res <- run_wilcoxon(sample_mp)
if (nrow(wilcox_res) > 0) {
  wilcox_res$feature_label <- ifelse(
    !is.na(mp_descriptions[wilcox_res$feature]),
    unname(mp_descriptions[wilcox_res$feature]),
    wilcox_res$feature
  )
}
write.csv(wilcox_res, "Auto_survival_tcga_cox_results.csv", row.names = FALSE)

plot_volcano <- function(df, title_text, out_pdf) {
  if (nrow(df) == 0) return(invisible(NULL))
  df <- df %>% mutate(sig = p_value < 0.05, neglog10 = -log10(p_value), effect_plot = effect)
  df$feature_plot <- ifelse(!is.na(mp_descriptions[df$feature]), unname(mp_descriptions[df$feature]), df$feature)

  p <- ggplot(df, aes(x = effect_plot, y = neglog10)) +
    geom_point(aes(color = sig), size = 3, alpha = 0.8) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick3"), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_text_repel(aes(label = feature_plot), size = 3.2, max.overlaps = 20) +
    theme_minimal(base_size = 13) +
    labs(title = title_text, x = "Median difference (Responder - Nonresponder)", y = "-log10(p)")

  ggsave(out_pdf, p, width = 9, height = 7)
}

plot_volcano(wilcox_res, "PDO response association volcano", "Auto_survival_tcga_volcano_EAC.pdf")

####################
# 4) Response-associated state analysis section
####################
if (length(state_groups) == 0) {
  state_groups <- list("All_MPs" = mp_cols)
}

state_names <- names(state_groups)
state_sample <- sample_mp
for (st in state_names) {
  use_mps <- intersect(state_groups[[st]], mp_cols)
  if (length(use_mps) == 0) next
  st_df <- sample_mp %>%
    filter(MP %in% use_mps) %>%
    group_by(orig.ident, Batch, Clinical.response.at.OG.MDT..responder.non.responder) %>%
    summarise(state_score = max(mp_mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(MP = st, mp_mean = state_score)
  state_sample <- bind_rows(state_sample, st_df %>% dplyr::select(colnames(sample_mp)))
}

state_only <- state_sample %>% filter(MP %in% state_names)
state_wilcox <- run_wilcoxon(state_only)
if (nrow(state_wilcox) > 0) {
  state_wilcox$feature_label <- state_wilcox$feature
}
write.csv(state_wilcox, "Auto_clinical_assoc_mps_states.csv", row.names = FALSE)

####################
# 5) Clinical association visualisation (single PDO cohort)
####################
assoc_df <- bind_rows(
  if (nrow(wilcox_res) > 0) wilcox_res %>% mutate(feature_type = "MP") else data.frame(),
  if (nrow(state_wilcox) > 0) state_wilcox %>% mutate(feature_type = "State") else data.frame()
)

plot_assoc <- function(df, out_pdf) {
  if (nrow(df) == 0) return(invisible(NULL))
  d <- df %>%
    mutate(
      neglog10p = -log10(p_value),
      sig = case_when(
        is.na(p_value) ~ "NS",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "NS"
      )
    )

  p <- ggplot(d, aes(x = feature, y = neglog10p)) +
    geom_col(aes(fill = effect), width = 0.75) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45", linewidth = 0.4) +
    geom_text(data = subset(d, sig != "NS"), aes(label = sig), color = "black", size = 3, vjust = -0.3) +
    facet_wrap(~ feature_type, scales = "free_x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Effect") +
    theme_minimal(base_size = 13) +
    labs(title = "PDO response-associated MPs and states", x = NULL, y = "-log10(p)") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 15),
      legend.position = "right"
    )

  ggsave(out_pdf, p, width = 14, height = 8)
}

plot_assoc(assoc_df, "Auto_clinical_assoc_EAC.pdf")

####################
# 6) Compact summary
####################
summary_df <- data.frame(
  metric = c(
    "n_pdo_samples",
    "n_mp_sets_tested",
    "n_response_assoc_mp_pval",
    "n_response_assoc_state_pval"
  ),
  value = c(
    dplyr::n_distinct(sample_mp$orig.ident),
    length(mp_cols),
    ifelse(nrow(wilcox_res) > 0, sum(wilcox_res$p_value < 0.05, na.rm = TRUE), 0),
    ifelse(nrow(state_wilcox) > 0, sum(state_wilcox$p_value < 0.05, na.rm = TRUE), 0)
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, "Auto_survival_clinical_mps_v2_summary.csv", row.names = FALSE)

cat("Saved: Auto_survival_tcga_volcano_EAC.pdf\n")
cat("Saved: Auto_survival_tcga_cox_results.csv\n")
cat("Saved: Auto_clinical_assoc_mps_states.csv\n")
cat("Saved: Auto_clinical_assoc_EAC.pdf\n")
cat("Saved: Auto_survival_clinical_mps_v2_summary.csv\n")
cat("=== Auto survival clinical MPs (PDO adapted) complete ===\n")
cat(format(Sys.time(), "%H:%M:%S"), "\n")
