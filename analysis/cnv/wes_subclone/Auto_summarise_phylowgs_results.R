####################
# Auto_summarise_phylowgs_results.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_summarise_phylowgs_results.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PhyloWGS *.summ.json.gz, *.muts.json.gz, and *.mutass.zip outputs.
#   - PhyloWGS input audit tables.
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/Auto_phylowgs_top_tree_summary.csv
#   - PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_phylowgs_*_top_tree.csv
# Downstream use: compact PhyloWGS tree/CNA assignment audit.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
live_root <- if (length(args) >= 1 && nzchar(args[1])) args[1] else Sys.getenv("OUT_ROOT")
if (!nzchar(live_root)) {
  live_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone"
}

samples <- c("PDO_1090_vs_NT_1090", "PDO_1181_vs_NT_1181")

read_json_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  fromJSON(txt = paste(readLines(con, warn = FALSE), collapse = "\n"), simplifyVector = FALSE)
}

as_num <- function(x) suppressWarnings(as.numeric(x))

children_for <- function(structure, pop_id) {
  kids <- structure[[as.character(pop_id)]]
  if (is.null(kids)) character()
  else as.character(unlist(kids, use.names = FALSE))
}

build_parent_map <- function(structure) {
  parent <- character()
  for (p in names(structure)) {
    kids <- as.character(unlist(structure[[p]], use.names = FALSE))
    if (length(kids)) parent[kids] <- p
  }
  parent
}

ancestor_chain <- function(pop_id, parent_map) {
  out <- as.character(pop_id)
  current <- as.character(pop_id)
  while (!is.na(parent_map[current]) && nzchar(parent_map[current])) {
    current <- parent_map[current]
    out <- c(current, out)
  }
  unique(out)
}

summary_rows <- list()

for (sample in samples) {
  sample_dir <- file.path(live_root, "reports/phylowgs", sample)
  table_dir <- file.path(live_root, "tables/phylowgs", sample)
  summ_path <- file.path(sample_dir, paste0(sample, ".summ.json.gz"))
  mutass_zip <- file.path(sample_dir, paste0(sample, ".mutass.zip"))
  cnv_audit_path <- file.path(table_dir, paste0("Auto_", sample, "_phylowgs_cnv_segment_audit.tsv"))
  ssm_map_path <- file.path(table_dir, paste0("Auto_", sample, "_phylowgs_ssm_pyclone_map.tsv"))
  if (!file.exists(summ_path)) stop("Missing PhyloWGS summary JSON: ", summ_path)
  if (!file.exists(mutass_zip)) stop("Missing PhyloWGS mutation-assignment zip: ", mutass_zip)
  if (!file.exists(cnv_audit_path)) stop("Missing CNV audit table: ", cnv_audit_path)
  if (!file.exists(ssm_map_path)) stop("Missing SSM/PyClone map: ", ssm_map_path)

  summ <- read_json_gz(summ_path)
  densities <- unlist(summ$tree_densities)
  densities <- densities[is.finite(as_num(densities))]
  if (!length(densities)) stop("No tree densities in ", summ_path)
  top_tree_id <- names(densities)[which.max(as_num(densities))]
  top_density <- as_num(densities[[top_tree_id]])
  top_density_fraction <- top_density / sum(as_num(densities))
  tree <- summ$trees[[top_tree_id]]
  if (is.null(tree)) stop("Top tree ID ", top_tree_id, " absent from ", summ_path)

  cnv_audit <- fread(cnv_audit_path)
  ssm_map <- fread(ssm_map_path)

  tmp_dir <- tempfile(paste0("Auto_", sample, "_mutass_"))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  mutass_file <- paste0(top_tree_id, ".json")
  unzip(mutass_zip, files = mutass_file, exdir = tmp_dir)
  mutass <- fromJSON(file.path(tmp_dir, mutass_file), simplifyVector = FALSE)
  assignments <- mutass$mut_assignments

  parent_map <- build_parent_map(tree$structure)
  pop_ids <- names(tree$populations)

  pop_tbl <- rbindlist(lapply(pop_ids, function(pop_id) {
    pop <- tree$populations[[pop_id]]
    data.table(
      sample = sample,
      top_tree_id = top_tree_id,
      population = pop_id,
      cellular_prevalence = paste(unlist(pop$cellular_prevalence), collapse = ";"),
      num_ssms = as.integer(pop$num_ssms),
      num_cnvs = as.integer(pop$num_cnvs),
      parent = ifelse(is.na(parent_map[pop_id]), NA_character_, parent_map[pop_id]),
      children = paste(children_for(tree$structure, pop_id), collapse = ";")
    )
  }), fill = TRUE)
  pop_tbl[, population_order := as.integer(population)]
  setorder(pop_tbl, population_order)
  pop_tbl[, population_order := NULL]

  cnv_assign <- rbindlist(lapply(names(assignments), function(pop_id) {
    cnvs <- as.character(unlist(assignments[[pop_id]]$cnvs, use.names = FALSE))
    if (!length(cnvs)) return(NULL)
    data.table(sample = sample, top_tree_id = top_tree_id, population = pop_id, cnv = cnvs)
  }), fill = TRUE)
  if (nrow(cnv_assign)) {
    cnv_assign <- merge(cnv_assign, cnv_audit, by = "cnv", all.x = TRUE)
    setcolorder(cnv_assign, intersect(c(
      "sample", "top_tree_id", "population", "cnv", "chrom", "start", "end",
      "major_cn", "minor_cn", "total_cn", "normal_cn", "cell_prev",
      "allele_cn_imputed", "n_overlapping_ssms"
    ), names(cnv_assign)))
    cnv_assign[, `:=`(population_order = as.integer(population), chrom_order = as.integer(chrom))]
    setorder(cnv_assign, population_order, chrom_order, start, end)
    cnv_assign[, `:=`(population_order = NULL, chrom_order = NULL)]
  }

  ssm_assign <- rbindlist(lapply(names(assignments), function(pop_id) {
    ssms <- as.character(unlist(assignments[[pop_id]]$ssms, use.names = FALSE))
    if (!length(ssms)) return(NULL)
    data.table(sample = sample, top_tree_id = top_tree_id, population = pop_id, ssm_id = ssms)
  }), fill = TRUE)
  if (nrow(ssm_assign)) {
    ssm_assign <- merge(ssm_assign, ssm_map, by = "ssm_id", all.x = TRUE)
    ssm_assign[, population_order := as.integer(population)]
    setorder(ssm_assign, population_order, chrom, pos)
    ssm_assign[, population_order := NULL]
  }

  clone_cna <- rbindlist(lapply(pop_ids, function(pop_id) {
    inherited_pops <- ancestor_chain(pop_id, parent_map)
    inherited_pops <- inherited_pops[inherited_pops != "0"]
    if (!nrow(cnv_assign) || !length(inherited_pops)) return(NULL)
    out <- copy(cnv_assign[population %in% inherited_pops])
    if (!nrow(out)) return(NULL)
    out[, `:=`(
      clone_population = pop_id,
      event_population = population,
      inherited_population_chain = paste(inherited_pops, collapse = ";")
    )]
    out
  }), fill = TRUE)
  if (nrow(clone_cna)) {
    setcolorder(clone_cna, intersect(c(
      "sample", "top_tree_id", "clone_population", "event_population",
      "inherited_population_chain", "cnv", "chrom", "start", "end",
      "major_cn", "minor_cn", "total_cn", "normal_cn", "cell_prev",
      "allele_cn_imputed", "n_overlapping_ssms"
    ), names(clone_cna)))
    clone_cna[, `:=`(clone_order = as.integer(clone_population), chrom_order = as.integer(chrom))]
    setorder(clone_cna, clone_order, chrom_order, start, end)
    clone_cna[, `:=`(clone_order = NULL, chrom_order = NULL)]
  }

  fwrite(pop_tbl, file.path(table_dir, paste0("Auto_", sample, "_phylowgs_population_top_tree.csv")))
  fwrite(cnv_assign, file.path(table_dir, paste0("Auto_", sample, "_phylowgs_cnv_assignment_top_tree.csv")))
  fwrite(ssm_assign, file.path(table_dir, paste0("Auto_", sample, "_phylowgs_ssm_assignment_top_tree.csv")))
  fwrite(clone_cna, file.path(table_dir, paste0("Auto_", sample, "_phylowgs_clone_cna_inherited_top_tree.csv")))

  summary_rows[[sample]] <- data.table(
    sample = sample,
    top_tree_id = top_tree_id,
    top_tree_density = top_density,
    top_tree_density_fraction = top_density_fraction,
    log_likelihood = as_num(tree$llh),
    n_posterior_trees = length(summ$trees),
    n_populations_including_root = length(pop_ids),
    n_event_populations = sum(pop_tbl$population != "0" & (pop_tbl$num_ssms > 0 | pop_tbl$num_cnvs > 0)),
    n_cnv_assigned_events = nrow(cnv_assign),
    n_ssm_assigned_events = nrow(ssm_assign),
    population_table = file.path(table_dir, paste0("Auto_", sample, "_phylowgs_population_top_tree.csv")),
    clone_cna_inherited_table = file.path(table_dir, paste0("Auto_", sample, "_phylowgs_clone_cna_inherited_top_tree.csv"))
  )
}

summary_tbl <- rbindlist(summary_rows, fill = TRUE)
out_path <- file.path(live_root, "tables/phylowgs", "Auto_phylowgs_top_tree_summary.csv")
fwrite(summary_tbl, out_path)
message("Wrote PhyloWGS top-tree summary: ", out_path)
