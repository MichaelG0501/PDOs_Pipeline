####################
# Analysis registry (authoritative):
#   Status: active terminal comparison; historical basename retained by request
#   Script: analysis/cell_states/legacy_Auto_pdo_sn_matched_pair_comparison.R
#   Methodology: analysis/methodology/cell_states/Auto_pdo_sn_matched_pair_comparison_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Compare matched PDO and biopsy current centred state/top-MP proportions.
#   Inputs under live PDOs_outs/:
#     PDOs_merged.rds; centred_mp_refinement/merged_refined_mp_genes.rds;
#     centred_mp_refinement/merged_refined_ucell_scores.rds;
#     centred_mp_refinement/centred_refined_noreg_{states,mp_adj}.rds;
#     centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv.
#   Inputs under live snSeq_Pipeline/sn_outs/:
#     snSeq_malignant_epi_meta.rds; cell_states/centred/intermediate/
#     merged_refined_mp_genes.rds; merged_refined_ucell_scores.rds;
#     centred_refined_noreg_{states,mp_adj}.rds;
#     cell_states/centred/tables/centred_refined_mp_state_grouping.csv.
#   Outputs under live PDOs_outs/Auto_pdo_sn_matched_pair_comparison/:
#     intermediate/Auto_matched_pair_cache.rds; tables/Auto_*.csv;
#     figures/Auto_final/01_current_centred_matched_pairs.pdf; reports/Auto_figure_index.csv;
#     logs/legacy_Auto_pdo_sn_matched_pair_comparison_run_summary.txt.
#   Downstream: none; descriptive terminal comparison for SUR680 and SUR791.
#   Cache: PDO_REPLOT_ONLY=1 requires unchanged input fingerprints; otherwise rebuild.
#   Run: qsub job_scripts/cell_states/legacy_Auto_pdo_sn_matched_pair_comparison.sh
#   Resources/environment: 2 CPUs, 32 GB, 02:00:00; dmtcp.
####################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_config.R")
source("/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/shared/Auto_pdo_analysis_helpers.R")
setwd(PDO_PROJECT_DIR)

####################
start_time <- Sys.time()
script_path <- "analysis/cell_states/legacy_Auto_pdo_sn_matched_pair_comparison.R"
out_dir <- Sys.getenv("PDO_SN_COMPARISON_OUT_DIR", file.path(PDO_LIVE_OUTS, "Auto_pdo_sn_matched_pair_comparison"))
paths <- pdo_ensure_output_tiers(out_dir)
cache_policy <- pdo_cache_policy()
cache_file <- file.path(paths["intermediate"], "Auto_matched_pair_cache.rds")
sn_root <- Sys.getenv("SNSEQ_OUTS_DIR", "/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/sn_outs")
sn_centred <- file.path(sn_root, "cell_states/centred")
pdo_centred <- file.path(PDO_LIVE_OUTS, "centred_mp_refinement")
input_files <- c(
  file.path(PDO_LIVE_OUTS, "PDOs_merged.rds"),
  file.path(pdo_centred, c("merged_refined_mp_genes.rds", "merged_refined_ucell_scores.rds",
                         "centred_refined_noreg_states.rds", "centred_refined_noreg_mp_adj.rds",
                         "tables/centred_refined_mp_state_grouping.csv")),
  file.path(sn_root, "snSeq_malignant_epi_meta.rds"),
  file.path(sn_centred, "intermediate", c("merged_refined_mp_genes.rds", "merged_refined_ucell_scores.rds",
                                        "centred_refined_noreg_states.rds", "centred_refined_noreg_mp_adj.rds")),
  file.path(sn_centred, "tables/centred_refined_mp_state_grouping.csv")
)
pdo_require_files(input_files)
# Size and mtime fingerprint avoids re-reading the 2.5-GB merged PDO file for replot.
info <- file.info(input_files)
fingerprints <- data.frame(path = input_files, size = info$size, mtime = as.numeric(info$mtime))
pair_map <- data.frame(patient_id = c("SUR680", "SUR791"),
                       pdo_sample = c("SUR680T3_PDO", "SUR791T3_PDO"),
                       sn_sample = c("H_post_T1_biopsy", "L_post_T1_biopsy"))

read_dataset <- function(meta, directory, intermediate, samples, modality, sample_column) {
  states <- pdo_normalise_final_state_vector(readRDS(file.path(intermediate, "centred_refined_noreg_states.rds")))
  adjusted <- as.matrix(readRDS(file.path(intermediate, "centred_refined_noreg_mp_adj.rds")))
  scores <- as.matrix(readRDS(file.path(intermediate, "merged_refined_ucell_scores.rds")))
  genes <- readRDS(file.path(intermediate, "merged_refined_mp_genes.rds"))
  grouping <- read.csv(file.path(directory, "tables/centred_refined_mp_state_grouping.csv"), check.names = FALSE)
  stopifnot("orig.ident" %in% names(meta), !is.null(names(states)), !anyDuplicated(names(states)),
            !anyDuplicated(rownames(meta)), !anyDuplicated(rownames(adjusted)),
            !anyDuplicated(rownames(scores)), !anyDuplicated(grouping$mp),
            all(c("mp", "state", "description") %in% names(grouping)),
            all(grouping$mp %in% names(genes)), all(grouping$mp %in% colnames(adjusted)),
            all(grouping$mp %in% colnames(scores)), !anyNA(meta$orig.ident))
  absent <- setdiff(samples, unique(as.character(meta$orig.ident)))
  if (length(absent)) stop(modality, " missing matched samples: ", paste(absent, collapse = ", "))
  cells <- rownames(meta)[meta$orig.ident %in% samples]
  if (!length(cells) || !all(cells %in% names(states)) || !all(cells %in% rownames(adjusted)) ||
      !all(cells %in% rownames(scores))) stop(modality, " incomplete state/score coverage for matched cells.")
  labels <- c(unique(grouping$state[grouping$state != "Cell cycle"]), "Unresolved", "Hybrid")
  if (!all(states[cells] %in% labels)) stop(modality, " has unexpected state labels.")
  adjusted <- adjusted[cells, grouping$mp, drop = FALSE]
  if (!all(is.finite(adjusted)) || !all(is.finite(scores[cells, grouping$mp, drop = FALSE]))) stop("Non-finite MP scores.")
  non_cc <- grouping$mp[grouping$state != "Cell cycle"]
  top_all <- grouping$mp[max.col(adjusted, ties.method = "first")]
  top_state <- non_cc[max.col(adjusted[, non_cc, drop = FALSE], ties.method = "first")]
  assignments <- data.frame(cell = cells, sample = as.character(meta[cells, "orig.ident"]),
                            modality = modality, state = unname(states[cells]),
                            top_mp = top_all, top_state_mp = top_state)
  assignments$patient_id <- pair_map$patient_id[match(assignments$sample, pair_map[[sample_column]])]
  assignments$top_state_mp_group <- grouping$state[match(top_state, grouping$mp)]
  state_props <- assignments %>% count(patient_id, sample, modality, state, name = "n") %>%
    complete(nesting(patient_id, sample, modality), state = labels, fill = list(n = 0L)) %>%
    group_by(patient_id, sample, modality) %>% mutate(total_n = sum(n), pct = 100 * n / total_n) %>% ungroup()
  mp_props <- bind_rows(lapply(c("top_mp", "top_state_mp"), function(column) {
    levels <- if (column == "top_mp") grouping$mp else non_cc
    assignments %>% mutate(mp = .data[[column]]) %>%
      count(patient_id, sample, modality, mp, name = "n") %>%
      complete(nesting(patient_id, sample, modality), mp = levels, fill = list(n = 0L)) %>%
      group_by(patient_id, sample, modality) %>% mutate(total_n = sum(n), pct = 100 * n / total_n) %>% ungroup() %>%
      mutate(selection = column) %>% left_join(grouping, by = "mp")
  }))
  list(assignments = assignments, states = state_props, mps = mp_props, grouping = grouping)
}

if (cache_policy$replot_only) {
  pdo_require_files(cache_file)
  cached <- readRDS(cache_file)
  if (!identical(cached$fingerprints, fingerprints)) stop("Inputs changed; rebuild before PDO_REPLOT_ONLY=1.")
  datasets <- cached$datasets
} else {
  pdo <- readRDS(file.path(PDO_LIVE_OUTS, "PDOs_merged.rds"))
  pdo_meta <- pdo[[]]
  rm(pdo)
  pdo_meta <- pdo_meta[as.character(pdo_meta$orig.ident) != PDO_EXCLUDED_SAMPLE, , drop = FALSE]
  datasets <- list(
    PDO = read_dataset(pdo_meta, pdo_centred, pdo_centred, pair_map$pdo_sample, "PDO organoid", "pdo_sample"),
    snRNA = read_dataset(readRDS(file.path(sn_root, "snSeq_malignant_epi_meta.rds")), sn_centred,
                        file.path(sn_centred, "intermediate"), pair_map$sn_sample, "snRNA-seq biopsy", "sn_sample")
  )
  saveRDS(list(datasets = datasets, fingerprints = fingerprints, pair_map = pair_map), cache_file)
}
state_props <- bind_rows(lapply(datasets, `[[`, "states"))
mp_props <- bind_rows(lapply(datasets, `[[`, "mps"))
assignments <- bind_rows(lapply(datasets, `[[`, "assignments"))
stopifnot(nrow(assignments) == sum(state_props$n))
write.csv(assignments, file.path(paths["tables"], "Auto_matched_pair_cell_assignments.csv"), row.names = FALSE)
write.csv(state_props, file.path(paths["tables"], "Auto_pdo_sn_matched_pair_state_proportions.csv"), row.names = FALSE)
write.csv(mp_props, file.path(paths["tables"], "Auto_pdo_sn_matched_pair_mp_proportions.csv"), row.names = FALSE)
summary <- state_props %>% group_by(patient_id, sample, modality) %>%
  summarise(cells = sum(n), dominant_state = state[which.max(n)], dominant_state_pct = max(pct), .groups = "drop")
write.csv(summary, file.path(paths["tables"], "Auto_pdo_sn_matched_pair_summary.csv"), row.names = FALSE)
colors <- c(PDO_STATE_COLORS, "Squamous-to-intestinal" = "#4DAF4A",
            "Glandular-to-intestinal" = "#FF7F00", "Cancer-cell immune mimicry" = "#377EB8", "Cell cycle" = "#777777")
# Separate native state labels and MP panels preserve dataset-specific meanings.
# Same MP identifiers in independent PDO and scRef collections are not equivalent.
final_dir <- file.path(paths["figures"], "Auto_final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(final_dir, "01_current_centred_matched_pairs.pdf")
pdf(pdf_path, width = 18, height = 12, onefile = TRUE)
for (patient in pair_map$patient_id) {
  panels <- lapply(names(datasets), function(dataset) {
    ds <- datasets[[dataset]]
    s <- ds$states[ds$states$patient_id == patient, ]
    m <- ds$mps[ds$mps$patient_id == patient & ds$mps$selection == "top_state_mp", ]
    s$state <- factor(s$state, levels = c(unique(ds$grouping$state[ds$grouping$state != "Cell cycle"]), "Unresolved", "Hybrid"))
    m$plot_label <- factor(m$plot_label, levels = ds$grouping$plot_label[ds$grouping$state != "Cell cycle"])
    p_state <- ggplot(s, aes(state, pct, fill = state)) + geom_col() + coord_flip() +
      scale_fill_manual(values = colors, guide = "none") +
      labs(x = NULL, y = "% of matched cells", title = paste(unique(s$sample), "native centred states")) + theme_bw(base_size = 16)
    p_mp <- ggplot(m, aes(plot_label, pct, fill = state)) + geom_col() + coord_flip() +
      scale_fill_manual(values = colors, guide = "none") +
      labs(x = NULL, y = "% of matched cells", title = "Top adjusted non-cell-cycle MP (all matched cells)") + theme_bw(base_size = 16)
    (p_state | p_mp) + plot_layout(widths = c(1, 1.5))
  })
  print(wrap_plots(panels, ncol = 1) +
          plot_annotation(title = paste(patient, "matched PDO and snRNA-seq biopsy"),
                          subtitle = "Independent native state definitions; MP IDs are dataset-specific. All cells remain in denominators."))
}
dev.off()
write.csv(data.frame(file = pdf_path, role = "Final descriptive comparison; one page per matched patient"),
          file.path(paths["reports"], "Auto_figure_index.csv"), row.names = FALSE)
pdo_write_run_summary(script_path, out_dir, inputs = input_files,
                      outputs = c(cache_file, pdf_path, list.files(paths["tables"], full.names = TRUE)),
                      parameters = list(start_time = as.character(start_time), end_time = as.character(Sys.time()),
                                        state_definitions = "current centred native; no 3CA relabel", patients = "SUR680,SUR791"),
                      cache = cache_policy)
message("MATCHED_PAIR_COMPARISON_OK: ", nrow(assignments), " matched cells")
####################
