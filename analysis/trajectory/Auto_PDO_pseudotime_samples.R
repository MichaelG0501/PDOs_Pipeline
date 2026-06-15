####################
# Auto_PDO_pseudotime_samples.R
#
# Per-sample Monocle3 pseudotime for PDOs using the four pre-relabel
# Approach-B states from Auto_PDO_states_noreg.rds. No batch correction is
# applied; each sample is preprocessed and ordered independently.
####################

helper_candidates <- c(
  "analysis/trajectory/Auto_PDO_pseudotime_helpers.R",
  "analysis/cell_states/Auto_PDO_pseudotime_helpers.R",
  "Auto_PDO_pseudotime_helpers.R"
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (is.na(helper_path)) stop("Could not find Auto_PDO_pseudotime_helpers.R")
source(helper_path)

message("=== PDO per-sample pre-relabel pseudotime ===")
results <- build_pdo_sample_trajectories(rebuild = TRUE)
if (length(results) == 0) stop("No sample trajectories were generated.")

pdf_path <- file.path(pseudotime_dir, "Auto_PDO_pseudotime_combined.pdf")
pdf(pdf_path, width = 14, height = 6, onefile = TRUE)
all_meta <- list()
for (sample_id in names(results)) {
  cds <- results[[sample_id]]$cds
  pt <- results[[sample_id]]$pseudotime
  meta_out <- results[[sample_id]]$metadata
  all_meta[[sample_id]] <- meta_out
  state_counts <- table(factor(meta_out$state, levels = primary_states))
  legend_labels <- setNames(paste0(primary_states, " (", as.integer(state_counts[primary_states]), ")"), primary_states)
  p_states <- plot_cells(
    cds,
    color_cells_by = "state_pre_relabel",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.8
  ) +
    scale_color_manual(values = state_cols[primary_states], breaks = primary_states, labels = legend_labels[primary_states], name = "State", na.value = "grey80", drop = FALSE, guide = guide_legend(override.aes = list(size = 4))) +
    labs(title = paste0("PDO states - ", sample_id, " (n = ", length(pt), ")"), color = NULL) +
    theme_minimal(base_size = 11)
  p_pseudotime <- plot_cells(
    cds,
    color_cells_by = "pseudotime",
    show_trajectory_graph = TRUE,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.8
  ) +
    scale_color_viridis_c(na.value = "grey85") +
    labs(title = paste0("Pseudotime - ", sample_id, " | root: ", root_state), color = "Pseudotime") +
    theme_minimal(base_size = 11)
  print(p_states + p_pseudotime + plot_layout(guides = "collect"))
}
dev.off()

meta_df <- bind_rows(all_meta)
summary_out <- meta_df %>%
  mutate(state = factor(state, levels = primary_states)) %>%
  group_by(sample, state) %>%
  summarise(n_cells = n(), n_cells_with_pseudotime = sum(is.finite(pseudotime)), median_pseudotime = safe_median(pseudotime), mean_pseudotime = safe_mean(pseudotime), .groups = "drop") %>%
  arrange(sample, state)
write_csv(meta_df, file.path(pseudotime_dir, "Auto_PDO_pseudotime_metadata.csv"))
write_csv(summary_out, file.path(pseudotime_dir, "Auto_PDO_pseudotime_summary.csv"))
write_csv(summary_out, file.path(summary_dir, "Auto_PDO_pseudotime_summary.csv"))
message("Saved PDO pseudotime outputs to: ", pseudotime_dir)
