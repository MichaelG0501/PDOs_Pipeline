####################
# Auto_PDO_numbat_concordance_summary_plots.R
#
# Presentation-style concordance summaries for InferCNA/arm-difference
# subclones versus Numbat haplotype-aware clone calls.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(scales)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

####################
# Match the heatmap workflow clone-mode switch. Conservative mode summarizes
# the re-cut robust Numbat clone layer without overwriting raw-summary plots.
clone_mode <- tolower(Sys.getenv("PDO_NUMBAT_CLONE_MODE", "raw"))
if (!clone_mode %in% c("raw", "conservative")) {
  stop("Unsupported PDO_NUMBAT_CLONE_MODE: ", clone_mode)
}
use_conservative_clones <- identical(clone_mode, "conservative")
numbat_label <- if (use_conservative_clones) "Numbat conservative" else "Numbat"
in_dir <- if (use_conservative_clones) "Auto_PDO_numbat/concordance_conservative" else "Auto_PDO_numbat/concordance"
out_pdf <- file.path(
  in_dir,
  if (use_conservative_clones) "Auto_PDO_numbat_conservative_concordance_summary_plots.pdf" else "Auto_PDO_numbat_concordance_summary_plots.pdf"
)
####################

required <- file.path(
  in_dir,
  c(
    "Auto_PDO_numbat_infercna_concordance_summary.csv",
    "Auto_PDO_numbat_infercna_contingency.csv",
    "Auto_PDO_numbat_clone_state_summary.csv",
    "Auto_PDO_numbat_clone_topmp_summary.csv"
  )
)
for (path in required) {
  if (!file.exists(path)) stop("Missing required concordance table: ", path)
}

summary_tbl <- fread(file.path(in_dir, "Auto_PDO_numbat_infercna_concordance_summary.csv"))
contingency_tbl <- fread(file.path(in_dir, "Auto_PDO_numbat_infercna_contingency.csv"))
state_tbl <- fread(file.path(in_dir, "Auto_PDO_numbat_clone_state_summary.csv"))
topmp_tbl <- fread(file.path(in_dir, "Auto_PDO_numbat_clone_topmp_summary.csv"))

state_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Basal to Intestinal Metaplasia" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved" = "grey80",
  "Hybrid" = "black"
)

mp_cols <- c(
  "MP6_G2M Cell Cycle" = "#B0B0B0",
  "MP7_DNA repair" = "#999999",
  "MP1_G2M checkpoint" = "#808080",
  "MP3_G1S Cell Cycle" = "#C0C0C0",
  "MP5_MYC-related Proliferation" = "#E41A1C",
  "MP4_Intestinal Metaplasia" = "#4DAF4A",
  "MP8_Columnar Progenitor" = "#FF7F00",
  "MP10_Inflammatory Stress Epi." = "#984EA3",
  "MP9_ECM Remodeling Epi." = "#C77CFF"
)

sample_order <- summary_tbl %>%
  mutate(rank_metric = ifelse(is.na(adjusted_rand), -Inf, adjusted_rand)) %>%
  arrange(desc(.data$rank_metric), desc(.data$normalised_mi), desc(.data$n_common_cells)) %>%
  pull(.data$sample)
summary_tbl$sample <- factor(summary_tbl$sample, levels = rev(sample_order))
contingency_tbl$sample <- factor(contingency_tbl$sample, levels = rev(sample_order))
state_tbl$sample <- factor(state_tbl$sample, levels = rev(sample_order))
topmp_tbl$sample <- factor(topmp_tbl$sample, levels = rev(sample_order))

metric_long <- summary_tbl %>%
  transmute(
    sample,
    ARI = adjusted_rand,
    NMI = normalised_mi,
    `Numbat purity within InferCNA` = purity_numbat_given_infercna,
    `InferCNA purity within Numbat` = purity_infercna_given_numbat
  ) %>%
  pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  mutate(metric = gsub("Numbat", numbat_label, .data$metric)) %>%
  mutate(value_plot = ifelse(is.na(.data$value), 0, .data$value),
         missing = is.na(.data$value))

p_metrics <- ggplot(metric_long, aes(.data$value_plot, .data$sample)) +
  geom_segment(aes(x = 0, xend = .data$value_plot, yend = .data$sample), color = "grey75", linewidth = 0.35) +
  geom_point(aes(alpha = !.data$missing), size = 2.1, color = "#2B6CB0") +
  facet_wrap(~metric, ncol = 2, scales = "free_x") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.25), guide = "none") +
  labs(title = paste0("InferCNA vs ", numbat_label, " clone concordance"), x = "Metric value", y = NULL) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey95", color = "grey80"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 7))

count_long <- summary_tbl %>%
  select(sample, n_infercna_subclones, n_numbat_clones) %>%
  pivot_longer(-sample, names_to = "method", values_to = "n_clones") %>%
  mutate(method = ifelse(.data$method == "n_numbat_clones", numbat_label, "InferCNA"))
p_counts <- ggplot(count_long, aes(.data$n_clones, .data$sample, color = .data$method)) +
  geom_line(aes(group = .data$sample), color = "grey80", linewidth = 0.3) +
  geom_point(size = 2.4) +
  scale_color_manual(values = setNames(c("#4C78A8", "#F58518"), c("InferCNA", numbat_label))) +
  labs(title = "Number of called subclones per sample", x = "Called subclones", y = NULL, color = NULL) +
  theme_classic(base_size = 10) +
  theme(axis.text.y = element_text(size = 7), legend.position = "top")

contingency_plot <- contingency_tbl %>%
  mutate(infercna_subclone = factor(.data$infercna_subclone, levels = sort(unique(.data$infercna_subclone))),
         numbat_clone = factor(.data$numbat_clone, levels = sort(unique(.data$numbat_clone))))
p_cont <- ggplot(contingency_plot, aes(.data$infercna_subclone, .data$numbat_clone, fill = .data$frac_of_infercna)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "#B2182B", labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  facet_wrap(~sample, scales = "free", ncol = 5) +
  labs(title = paste0(numbat_label, " clone composition inside each InferCNA subclone"), x = "InferCNA subclone", y = paste0(numbat_label, " clone"), fill = "% of InferCNA") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 5.5))

state_plot <- state_tbl %>%
  mutate(state_label = ifelse(is.na(.data$state_label), "NA", as.character(.data$state_label)),
         state_label = factor(.data$state_label, levels = c(names(state_cols), "NA")))
state_cols_use <- state_cols[names(state_cols) %in% levels(state_plot$state_label)]
state_cols_use <- c(state_cols_use, "NA" = "grey90")
p_state <- ggplot(state_plot, aes(.data$numbat_clone, .data$frac, fill = .data$state_label)) +
  geom_col(color = "grey20", linewidth = 0.08) +
  scale_fill_manual(values = state_cols_use, drop = FALSE, labels = function(x) gsub("Basal to Intestinal Metaplasia", "Basal to Intest. Meta", x)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  facet_wrap(~sample, scales = "free_x", ncol = 5) +
  labs(title = paste0("Final PDO state composition of ", numbat_label, " clones"), x = paste0(numbat_label, " clone"), y = "% cells", fill = "State") +
  theme_classic(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 6),
        legend.position = "bottom")

topmp_plot <- topmp_tbl %>%
  mutate(top_mp_label = ifelse(is.na(.data$top_mp_label), "NA", as.character(.data$top_mp_label)))
mp_cols_use <- mp_cols[names(mp_cols) %in% unique(topmp_plot$top_mp_label)]
missing_mp <- setdiff(unique(topmp_plot$top_mp_label), names(mp_cols_use))
if (length(missing_mp) > 0) {
  mp_cols_use <- c(mp_cols_use, setNames(hue_pal()(length(missing_mp)), missing_mp))
}
p_topmp <- ggplot(topmp_plot, aes(.data$numbat_clone, .data$frac, fill = .data$top_mp_label)) +
  geom_col(color = "grey20", linewidth = 0.08) +
  scale_fill_manual(values = mp_cols_use, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  facet_wrap(~sample, scales = "free_x", ncol = 5) +
  labs(title = paste0("Top PDO metaprogram composition of ", numbat_label, " clones"), x = paste0(numbat_label, " clone"), y = "% cells", fill = "Top MP") +
  theme_classic(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        axis.text.y = element_text(size = 6),
        legend.position = "bottom")

pdf(out_pdf, width = 15, height = 10, useDingbats = FALSE)
grid.arrange(p_metrics, p_counts, ncol = 2, widths = c(1.25, 1))
print(p_cont)
print(p_state)
print(p_topmp)
dev.off()

message("Wrote: ", file.path(out_root, out_pdf))
