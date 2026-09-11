####################
# Analysis registry:
#   Status: active terminal clinical and FLOT-sensitivity association workflow
#   Script: analysis/clinical/Auto_centred_clinical_auc_maxstat.R
#   Methodology: none
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Finds the optimal cutpoint of AUCrel that maximizes the difference
#     in MP UCell score or state proportion using the maxstat package.
#     Also plots boxplots of the optimal splits and bubble plots of all possible splits.
#   Inputs:
#     PDOs_outs/Auto_centred_clinical_auc_response_associations/intermediate/Auto_centred_clinical_auc_response_results.rds
#   Outputs:
#     PDOs_outs/Auto_centred_clinical_auc_maxstat/tables/*
#     PDOs_outs/Auto_centred_clinical_auc_maxstat/figures/*
#   Downstream: none
#   Cache/replot: none
#   Run: qsub analysis/clinical/Auto_centred_clinical_auc_maxstat.sh
#   Environment: /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
####################

suppressPackageStartupMessages({
  library(tidyverse)
  library(maxstat)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

project_dir <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_config.R"))
source(file.path(project_dir, "analysis/shared/Auto_pdo_analysis_helpers.R"))

in_dir <- file.path(PDO_OUTPUT_DIR, "Auto_centred_clinical_auc_response_associations")
cache_path <- file.path(in_dir, "intermediate", "Auto_centred_clinical_auc_response_results.rds")

out_dir <- file.path(PDO_OUTPUT_DIR, "Auto_centred_clinical_auc_maxstat")
out_paths <- pdo_ensure_output_tiers(out_dir)
tables_dir <- out_paths[["tables"]]
figures_dir <- out_paths[["figures"]]

if (!file.exists(cache_path)) {
  stop("Input cache not found: ", cache_path)
}

res <- readRDS(cache_path)
auc_state <- res$auc_state
auc_mp <- res$auc_mp

state_order <- res$state_order
mp_order <- res$mp_order

state_axis_labels <- c(
  "Classic proliferation" = "Classic\nproliferation",
  "Columnar-to-intestinal" = "Columnar-to-\nintestinal",
  "Glandular differentiation" = "Glandular\ndifferentiation",
  "Stress-adaptive" = "Stress-\nadaptive",
  "ECM-remodelling" = "ECM-\nremodelling",
  "Motile-cilia differentiation" = "Motile-cilia\ndifferentiation"
)

mp_labels <- setNames(
  ifelse(is.na(PDO_MP_DESCRIPTIONS[mp_order]), mp_order, paste0(mp_order, "\n", PDO_MP_DESCRIPTIONS[mp_order])),
  mp_order
)


run_maxstat_and_all_cuts <- function(data, feature_col, is_state=FALSE) {
  features <- unique(data[[feature_col]])
  
  opt_results <- list()
  all_cuts_results <- list()
  
  for (f in features) {
    df <- data[data[[feature_col]] == f, ]
    df <- df[!is.na(df$auc_rel) & !is.na(df$value), ]
    
    if (nrow(df) < 6) next
    
    # 1. Optimal cut using maxstat
    test_res <- maxstat.test(value ~ auc_rel, data = df, 
                             smethod = "Wilcoxon", pmethod = "exact", 
                             minprop = 0.25, maxprop = 0.75)
                             
    cutpoint <- test_res$estimate
    pval <- test_res$p.value
    
    high_group <- df$value[df$auc_rel > cutpoint]
    low_group <- df$value[df$auc_rel <= cutpoint]
    n_high <- length(high_group)
    n_low <- length(low_group)
    
    eff_size <- median(high_group, na.rm=TRUE) - median(low_group, na.rm=TRUE)
    
    opt_results[[f]] <- data.frame(
      feature = f,
      optimal_cutpoint = cutpoint,
      p_value = pval,
      effect_size = eff_size,
      n_high = n_high,
      n_low = n_low,
      stringsAsFactors = FALSE
    )
    
    # 2. All possible cuts within 0.25 and 0.75 quantiles
    sorted_auc <- sort(unique(df$auc_rel))
    n <- length(sorted_auc)
    # the cuts are typically the sorted values (or midpoints). 
    # To get exactly the splits Top 3 to Top 9, we take the 3rd to 9th largest AUCrel.
    # largest AUCrel is sorted_auc[12]. 3rd largest is sorted_auc[10].
    # So we loop over possible split points that produce 3 to 9 samples in high_group.
    
    possible_cuts <- sorted_auc[3:9] # Wait, if n=12, sorted_auc[3] means 9 samples are > cutpoint.
    
    cuts_list <- lapply(sorted_auc, function(c) {
      g_high <- df$value[df$auc_rel > c]
      g_low <- df$value[df$auc_rel <= c]
      
      nh <- length(g_high)
      nl <- length(g_low)
      if (nh < 3 || nh > 9) return(NULL) # keep min 3 per group
      
      wt <- suppressWarnings(wilcox.test(g_high, g_low, exact=FALSE))
      e_size <- median(g_high, na.rm=TRUE) - median(g_low, na.rm=TRUE)
      
      split_name <- paste0("Top ", nh, "\nvs\nBottom ", nl)
      
      data.frame(
        feature = f,
        cutpoint = c,
        split_name = split_name,
        n_high = nh,
        p_value = wt$p.value,
        effect_size = e_size,
        stringsAsFactors = FALSE
      )
    })
    cuts_list <- Filter(Negate(is.null), cuts_list)
    if(length(cuts_list) > 0) all_cuts_results[[f]] <- do.call(rbind, cuts_list)
  }
  
  opt_df <- do.call(rbind, opt_results)
  if(is_state) {
    opt_df$feature <- factor(opt_df$feature, levels = state_order)
  } else {
    opt_df$feature <- factor(opt_df$feature, levels = mp_order)
  }
  
  list(
    optimal = opt_df,
    all_cuts = do.call(rbind, all_cuts_results)
  )
}

state_res <- run_maxstat_and_all_cuts(auc_state, "feature", is_state=TRUE)
mp_res <- run_maxstat_and_all_cuts(auc_mp, "feature", is_state=FALSE)

state_optimal <- state_res$optimal
mp_optimal <- mp_res$optimal

state_all <- state_res$all_cuts
mp_all <- mp_res$all_cuts

# Fix split_name ordering for bubble plot
split_levels <- paste0("Top ", 9:3, "\nvs\nBottom ", 3:9)
state_all$split_name <- factor(state_all$split_name, levels = split_levels)
mp_all$split_name <- factor(mp_all$split_name, levels = split_levels)

# Map labels to features
state_optimal$feature_label <- state_axis_labels[as.character(state_optimal$feature)]
mp_optimal$feature_label <- mp_labels[as.character(mp_optimal$feature)]
state_all$feature_label <- state_axis_labels[as.character(state_all$feature)]
mp_all$feature_label <- mp_labels[as.character(mp_all$feature)]

state_all$feature_label <- factor(state_all$feature_label, levels = state_axis_labels[state_order])
mp_all$feature_label <- factor(mp_all$feature_label, levels = mp_labels[mp_order])

# adjust p-values for optimal
state_optimal$padj <- p.adjust(state_optimal$p_value, method="BH")
mp_optimal$padj <- p.adjust(mp_optimal$p_value, method="BH")

write.csv(state_optimal, file.path(tables_dir, "Auto_maxstat_state_optimal.csv"), row.names=FALSE)
write.csv(mp_optimal, file.path(tables_dir, "Auto_maxstat_mp_optimal.csv"), row.names=FALSE)

# ---------------------------------------------------------
# Plot 1: Volcano plot
# ---------------------------------------------------------
plot_volcano <- function(res, title, y_label, x_label="Effect Size (Median Diff: High - Low AUCrel)") {
  res$logP <- -log10(res$p_value)
  res$sig <- res$p_value < 0.05
  
  ggplot(res, aes(x = effect_size, y = logP)) +
    geom_point(aes(color = sig), size=3, alpha=0.8) +
    geom_text_repel(aes(label = feature_label), size=2.5, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.05), linetype="dashed", color="red") +
    geom_vline(xintercept = 0, linetype="dashed", color="black") +
    scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "red")) +
    theme_classic() +
    labs(title = title, x = x_label, y = y_label, color = "P < 0.05")
}

p_state_volc <- plot_volcano(state_optimal, "Maxstat Optimal Cutpoint (State Proportions)", "-log10(P-value)")
p_mp_volc <- plot_volcano(mp_optimal, "Maxstat Optimal Cutpoint (MP UCell Scores)", "-log10(P-value)")

ggsave(file.path(figures_dir, "Auto_maxstat_state_volcano.pdf"), p_state_volc, width = 9, height = 7)
ggsave(file.path(figures_dir, "Auto_maxstat_mp_volcano.pdf"), p_mp_volc, width = 11, height = 8)

# ---------------------------------------------------------
# Plot 2: Boxplots for the optimal cuts in a single row
# ---------------------------------------------------------
plot_optimal_boxplots_single_row <- function(data, opt_res, label_map, title, y_label, percent_axis = FALSE) {
  # Merge data with optimal cutpoint
  plot_data <- data %>%
    inner_join(opt_res %>% select(feature, optimal_cutpoint, p_value, n_high, n_low), by="feature") %>%
    mutate(
      Group = ifelse(auc_rel > optimal_cutpoint, "High AUCrel (Resistant)", "Low AUCrel (Sensitive)"),
      Group = factor(Group, levels = c("Low AUCrel (Sensitive)", "High AUCrel (Resistant)"))
    )
  
  # Ensure feature ordering
  plot_data$feature <- factor(plot_data$feature, levels = levels(opt_res$feature))
  
  # Format x-axis labels with split count
  # e.g., "Feature Name\n(Top 3 vs 9)"
  custom_labels <- sapply(levels(opt_res$feature), function(f) {
    row <- opt_res[opt_res$feature == f, ]
    if(nrow(row) == 0) return(f)
    # label_map maps feature string to pretty name
    base_name <- label_map[f]
    paste0(base_name, "\n(Top ", row$n_high, " vs Bottom ", row$n_low, ")")
  })
  
  # P-value annotations
  annotation <- opt_res %>%
    group_by(feature) %>%
    summarise(
      y_max = max(data$value[data$feature == feature], na.rm = TRUE),
      y_min = min(data$value[data$feature == feature], na.rm = TRUE),
      p_value = p_value[1],
      .groups = "drop"
    ) %>%
    mutate(
      y = y_max + pmax(ifelse(percent_axis, 2, 0.01), 0.1 * pmax(y_max - y_min, ifelse(percent_axis, 5, 0.02))),
      label = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
    
  p <- ggplot(plot_data, aes(x = feature, y = value, fill = Group)) +
    geom_boxplot(position = position_dodge(width = 0.75), width = 0.62, outlier.shape = NA, alpha = 0.8, colour = "black") +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75), size = 1.5, shape = 16, colour = "#222222", alpha=0.9, stroke = 0) +
    geom_text(data = annotation, aes(x = feature, y = y, label = label), inherit.aes = FALSE, size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c("Low AUCrel (Sensitive)" = "#0072B2", "High AUCrel (Resistant)" = "#D55E00")) +
    scale_x_discrete(labels = custom_labels) +
    labs(title = title, x = NULL, y = y_label, fill = "Optimal AUCrel Split") +
    pdo_theme_slide(11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
    coord_cartesian(clip = "off")
    
  if (percent_axis) p <- p + scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0.02, 0.18)))
  else p <- p + scale_y_continuous(expand = expansion(mult = c(0.02, 0.18)))
  p
}

p_state_box <- plot_optimal_boxplots_single_row(auc_state, state_optimal, state_axis_labels, "Optimal maxstat Split: State Proportions", "Sample state proportion (%)", TRUE)
p_mp_box <- plot_optimal_boxplots_single_row(auc_mp, mp_optimal, mp_labels, "Optimal maxstat Split: MP UCell Scores", "Mean sample UCell score", FALSE)

ggsave(file.path(figures_dir, "Auto_maxstat_state_optimal_boxplots.pdf"), p_state_box, width = 14, height = 7.5)
ggsave(file.path(figures_dir, "Auto_maxstat_mp_optimal_boxplots.pdf"), p_mp_box, width = 18, height = 8.5)


# ---------------------------------------------------------
# Plot 3: Bubble plots for ALL evaluated cutpoints
# ---------------------------------------------------------
plot_bubble <- function(all_cuts_data, title) {
  all_cuts_data$logP <- -log10(all_cuts_data$p_value)
  all_cuts_data$logP[all_cuts_data$logP > 4] <- 4
  
  ggplot(all_cuts_data, aes(x = split_name, y = feature_label)) +
    geom_point(aes(size = logP, fill = effect_size), shape=21, color="black", alpha=0.85) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, name = "Effect Size\n(High - Low)") +
    scale_size_continuous(range = c(2, 8), name = "-log10(P-value)\n(Wilcoxon)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          panel.grid.major.x = element_line(color="grey80", linetype="dashed")) +
    labs(title = title, x = "Evaluated Split (High AUCrel vs Low AUCrel)", y = "Feature")
}

p_state_bub <- plot_bubble(state_all, "Landscape of Splits (State Proportions)")
p_mp_bub <- plot_bubble(mp_all, "Landscape of Splits (MP UCell Scores)")

ggsave(file.path(figures_dir, "Auto_maxstat_state_all_cuts_bubble.pdf"), p_state_bub, width = 10, height = 7)
ggsave(file.path(figures_dir, "Auto_maxstat_mp_all_cuts_bubble.pdf"), p_mp_bub, width = 12, height = 10)

message("Maxstat script completed successfully with all visualisations.")
