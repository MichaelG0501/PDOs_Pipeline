####################
# Auto_five_state_surface_markers.R
#
# Prioritize FACS-suitable surface markers for the five finalized PDO states.
#
# NOTE:
#   This script downloads and caches two external annotation references if they
#   are not already present:
#     1) UniProt reviewed human subcellular-location + topology table
#     2) ETH Zurich human surfaceome Table S3 workbook
#   Default cache directory:
#     /rds/general/project/spatialtranscriptomics/ephemeral/Auto_pdo_surface_marker_db
#   Fallback cache directory:
#     PDOs_outs/Auto_five_state_surface_markers/db_cache
#
# Inputs:
#   PDOs_outs/PDOs_merged.rds
#   PDOs_outs/Auto_PDO_final_states.rds
#   PDOs_outs/Auto_five_state_markers/Auto_five_state_marker_summary.csv
#   PDOs_outs/Auto_five_state_markers/Auto_five_state_markers_ranked.csv
#
# Outputs:
#   PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_ranked.csv
#   PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_state_metrics.csv
#   PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_database_manifest.csv
#   PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_candidates.xlsx
####################

####################
# libraries
####################
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(openxlsx)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

####################
# setup
####################
project_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
setwd(file.path(project_dir, "PDOs_outs"))

marker_dir <- "Auto_five_state_markers"
out_dir <- "Auto_five_state_surface_markers"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

state_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)

state_cols <- c(
  "Classic Proliferative" = "#E41A1C", # Red
  "Basal to Intest. Meta" = "#4DAF4A", # Green
  "SMG-like Metaplasia"   = "#FF7F00", # Orange
  "Stress-adaptive"       = "#984EA3", # Purple
  "3CA_EMT_and_Protein_maturation" = "#377EB8"  # Blue
)

sheet_map <- c(
  "Classic Proliferative" = "Classic_Prolif",
  "Basal to Intest. Meta" = "Basal_Int_Meta",
  "Stress-adaptive" = "Stress_adapt",
  "SMG-like Metaplasia" = "SMG_metaplasia",
  "3CA_EMT_and_Protein_maturation" = "3CA_EMT_ProtMat"
)

uniprot_url <- paste0(
  "https://rest.uniprot.org/uniprotkb/stream?format=tsv&query=",
  utils::URLencode("(organism_id:9606) AND (reviewed:true)", reserved = TRUE),
  "&fields=accession,gene_primary,protein_name,cc_subcellular_location,ft_transmem,ft_topo_dom"
)

surfaceome_url <- "https://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx"

####################
# helpers
####################
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  path
}

resolve_db_dir <- function(primary_dir, fallback_dir) {
  ok <- tryCatch({
    ensure_dir(primary_dir)
    file.access(primary_dir, mode = 2) == 0
  }, error = function(e) FALSE)

  if (isTRUE(ok)) {
    return(primary_dir)
  }

  ensure_dir(fallback_dir)
}

download_if_missing <- function(url, dest) {
  ensure_dir(dirname(dest))

  downloaded_now <- FALSE
  if (!file.exists(dest)) {
    message("Downloading annotation reference: ", basename(dest))
    utils::download.file(url = url, destfile = dest, mode = "wb", quiet = FALSE)
    downloaded_now <- TRUE
  }

  list(
    path = dest,
    downloaded_now = downloaded_now,
    file_size = file.info(dest)$size %||% NA_real_,
    file_mtime = as.character(file.info(dest)$mtime %||% NA)
  )
}

collapse_unique <- function(x, sep = " || ") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

expand_go_terms <- function(root_ids) {
  offspring <- unlist(lapply(root_ids, function(go_id) {
    if (exists(go_id, envir = GOCCOFFSPRING, inherits = FALSE)) {
      get(go_id, envir = GOCCOFFSPRING)
    } else {
      NA_character_
    }
  }))

  out <- unique(c(root_ids, offspring))
  out[!is.na(out)]
}

get_assay_layer <- function(obj, assay = "RNA", layer = "data") {
  out <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- tryCatch(
      GetAssayData(obj, assay = assay, slot = layer),
      error = function(e) NULL
    )
  }

  if (is.null(out)) {
    stop("Could not extract assay layer/slot '", layer, "'.")
  }

  out
}

safe_sheet_name <- function(x) {
  out <- sheet_map[[x]] %||% gsub("[^A-Za-z0-9_]", "_", x)
  substr(out, 1, 31)
}

make_flag_summary <- function(df) {
  apply(
    df,
    1,
    function(row_vals) {
      labels <- names(row_vals)[as.logical(row_vals)]
      if (length(labels) == 0) "none" else paste(labels, collapse = "; ")
    }
  )
}

percent_rank0 <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  dplyr::percent_rank(dplyr::coalesce(x, min(x, na.rm = TRUE)))
}

####################
# database cache setup
####################
db_dir <- resolve_db_dir(
  primary_dir = Sys.getenv(
    "PDO_SURFACE_DB_DIR",
    unset = "/rds/general/project/spatialtranscriptomics/ephemeral/Auto_pdo_surface_marker_db"
  ),
  fallback_dir = file.path(getwd(), out_dir, "db_cache")
)

message("Using annotation cache directory: ", db_dir)

uniprot_ref <- download_if_missing(
  url = uniprot_url,
  dest = file.path(db_dir, "Auto_uniprot_human_reviewed_surface_topology.tsv")
)

surfaceome_ref <- download_if_missing(
  url = surfaceome_url,
  dest = file.path(db_dir, "Auto_table_S3_surfaceome.xlsx")
)

db_manifest <- data.frame(
  source = c(
    "UniProt reviewed human subcellular location + topology",
    "ETH Zurich human surfaceome Table S3"
  ),
  url = c(uniprot_url, surfaceome_url),
  local_file = c(uniprot_ref$path, surfaceome_ref$path),
  downloaded_now = c(uniprot_ref$downloaded_now, surfaceome_ref$downloaded_now),
  file_size_bytes = c(uniprot_ref$file_size, surfaceome_ref$file_size),
  file_mtime = c(uniprot_ref$file_mtime, surfaceome_ref$file_mtime),
  stringsAsFactors = FALSE
)

fwrite(
  db_manifest,
  file.path(out_dir, "Auto_five_state_surface_marker_database_manifest.csv")
)

####################
# marker inputs
####################
message("Loading five-state PDO marker summaries.")

marker_summary <- fread(file.path(marker_dir, "Auto_five_state_marker_summary.csv")) %>%
  as_tibble()

ranked_markers <- fread(file.path(marker_dir, "Auto_five_state_markers_ranked.csv")) %>%
  as_tibble() %>%
  dplyr::select(
    state,
    gene,
    reproducibility_rank,
    effect_rank,
    specificity_rank,
    ranking_score
  )

candidate_tbl <- marker_summary %>%
  left_join(ranked_markers, by = c("state", "gene")) %>%
  filter(
    hit_sample_n > 0,
    best_state_match,
    specificity_gap > 0
  ) %>%
  distinct(state, gene, .keep_all = TRUE)

if (nrow(candidate_tbl) == 0) {
  stop("No five-state PDO markers passed the initial reuse filter.")
}

####################
# current UniProt annotation
####################
message("Parsing UniProt reviewed-human annotations.")

uniprot_raw <- fread(
  uniprot_ref$path,
  sep = "\t",
  header = TRUE,
  quote = ""
) %>%
  as_tibble() %>%
  transmute(
    gene = trimws(`Gene Names (primary)`),
    uniprot_entry = Entry,
    uniprot_protein = `Protein names`,
    uniprot_subcellular = `Subcellular location [CC]`,
    uniprot_transmem = Transmembrane,
    uniprot_topology = `Topological domain`
  ) %>%
  filter(!is.na(gene), nzchar(gene))

surface_regex <- paste(
  c(
    "cell membrane",
    "plasma membrane",
    "cell surface",
    "apical cell membrane",
    "basolateral cell membrane",
    "lateral cell membrane",
    "basal cell membrane"
  ),
  collapse = "|"
)

uniprot_gene <- uniprot_raw %>%
  group_by(gene) %>%
  summarise(
    uniprot_entry = collapse_unique(uniprot_entry),
    uniprot_protein = collapse_unique(uniprot_protein),
    uniprot_subcellular = collapse_unique(uniprot_subcellular),
    uniprot_transmem = collapse_unique(uniprot_transmem),
    uniprot_topology = collapse_unique(uniprot_topology),
    uniprot_surface_flag = grepl(
      surface_regex,
      tolower(uniprot_subcellular %||% ""),
      perl = TRUE
    ),
    uniprot_extracellular_flag = grepl(
      "extracellular",
      uniprot_topology %||% "",
      ignore.case = TRUE
    ),
    uniprot_transmem_flag = grepl(
      "TRANSMEM",
      uniprot_transmem %||% "",
      fixed = TRUE
    ),
    .groups = "drop"
  )

####################
# surfaceome reference
####################
message("Parsing surfaceome workbook.")

surfaceome_raw <- suppressWarnings(
  readxl::read_excel(surfaceome_ref$path, skip = 1)
) %>%
  as_tibble() %>%
  transmute(
    gene = trimws(`UniProt gene`),
    surfaceome_label = `Surfaceome Label`,
    surfaceome_label_source = `Surfaceome Label Source`,
    surfaceome_tm_domains = suppressWarnings(as.numeric(`TM domains`)),
    surfaceome_topology = topology,
    cspa_category = `CSPA category`,
    surfaceome_uniprot_subcellular = `UniProt subcellular`
  ) %>%
  filter(!is.na(gene), nzchar(gene))

surfaceome_gene <- surfaceome_raw %>%
  group_by(gene) %>%
  summarise(
    surfaceome_label = collapse_unique(surfaceome_label),
    surfaceome_label_source = collapse_unique(surfaceome_label_source),
    surfaceome_topology = collapse_unique(surfaceome_topology),
    cspa_category = collapse_unique(cspa_category),
    surfaceome_uniprot_subcellular = collapse_unique(surfaceome_uniprot_subcellular),
    surfaceome_surface_flag = any(tolower(surfaceome_label) == "surface", na.rm = TRUE),
    surfaceome_extracellular_flag = any(grepl("NC:", surfaceome_topology, fixed = TRUE), na.rm = TRUE),
    surfaceome_transmem_flag = any(surfaceome_tm_domains > 0, na.rm = TRUE),
    cspa_support_flag = any(!is.na(cspa_category) & nzchar(cspa_category), na.rm = TRUE),
    cspa_high_confidence_flag = any(grepl("high confidence", cspa_category, ignore.case = TRUE), na.rm = TRUE),
    .groups = "drop"
  )

####################
# GO cellular-component support
####################
message("Building GO plasma-membrane/cell-surface support.")

go_surface_terms <- expand_go_terms(c("GO:0009986", "GO:0009897"))
go_plasma_terms <- expand_go_terms(c("GO:0005886", "GO:0005887"))
go_external_terms <- expand_go_terms(c("GO:0009897"))

candidate_genes <- sort(unique(candidate_tbl$gene))

go_raw <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = candidate_genes,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GOALL", "ONTOLOGYALL")
) %>%
  as_tibble()

go_gene <- go_raw %>%
  filter(ONTOLOGYALL == "CC", !is.na(GOALL)) %>%
  group_by(SYMBOL) %>%
  summarise(
    go_cc_terms = collapse_unique(GOALL, sep = "; "),
    go_surface_flag = any(GOALL %in% go_surface_terms),
    go_plasma_membrane_flag = any(GOALL %in% go_plasma_terms),
    go_external_side_flag = any(GOALL %in% go_external_terms),
    .groups = "drop"
  ) %>%
  dplyr::rename(gene = SYMBOL)

####################
# expression summaries
####################
message("Computing sample-level median expression summaries from PDOs_merged.rds.")

pdos_all <- readRDS("PDOs_merged.rds")
state_labels <- readRDS("Auto_PDO_final_states.rds")

DefaultAssay(pdos_all) <- "RNA"

common_cells <- intersect(colnames(pdos_all), names(state_labels))
state_labels <- state_labels[common_cells]
keep_cells <- common_cells[as.character(state_labels) %in% state_order]

candidate_genes <- intersect(candidate_genes, rownames(pdos_all))
candidate_tbl <- candidate_tbl %>%
  filter(gene %in% candidate_genes)

if (length(candidate_genes) == 0) {
  stop("No candidate marker genes were found in PDOs_merged.rds.")
}

counts_mat <- get_assay_layer(pdos_all, assay = "RNA", layer = "counts")[candidate_genes, keep_cells, drop = FALSE]
data_mat <- get_assay_layer(pdos_all, assay = "RNA", layer = "data")[candidate_genes, keep_cells, drop = FALSE]

sample_labels <- as.character(pdos_all@meta.data$orig.ident[match(keep_cells, colnames(pdos_all))])
unique_samples <- unique(sample_labels)

state_factor <- factor(as.character(state_labels[keep_cells]), levels = state_order)
state_cells <- split(keep_cells, state_factor)

pct_mat <- matrix(NA_real_, nrow = length(candidate_genes), ncol = length(state_order))
mean_mat <- matrix(NA_real_, nrow = length(candidate_genes), ncol = length(state_order))
rownames(pct_mat) <- candidate_genes
colnames(pct_mat) <- state_order
rownames(mean_mat) <- candidate_genes
colnames(mean_mat) <- state_order

for (state_name in state_order) {
  cells_in_state <- state_cells[[state_name]] %||% character()
  if (length(cells_in_state) == 0) next
  
  sample_pcts <- list()
  sample_means <- list()
  
  for (samp in unique_samples) {
    cells_samp <- cells_in_state[sample_labels[match(cells_in_state, keep_cells)] == samp]
    if (length(cells_samp) > 0) {
      samp_counts <- counts_mat[, cells_samp, drop = FALSE]
      samp_data <- data_mat[, cells_samp, drop = FALSE]
      
      sample_pcts[[samp]] <- Matrix::rowMeans(samp_counts > 0)
      sample_means[[samp]] <- Matrix::rowMeans(samp_data)
    }
  }
  
  if (length(sample_pcts) > 0) {
    pct_mat[, state_name] <- apply(do.call(cbind, sample_pcts), 1, median, na.rm = TRUE)
    mean_mat[, state_name] <- apply(do.call(cbind, sample_means), 1, median, na.rm = TRUE)
  }
}

state_metric_long <- bind_rows(lapply(state_order, function(state_name) {
  data.frame(
    gene = candidate_genes,
    state = state_name,
    pct_cells = pct_mat[, state_name],
    mean_logexpr = mean_mat[, state_name],
    stringsAsFactors = FALSE
  )
}))

fwrite(
  state_metric_long,
  file.path(out_dir, "Auto_five_state_surface_marker_state_metrics.csv")
)

####################
# merge annotation + expression
####################
message("Scoring FACS suitability.")

scored_tbl <- candidate_tbl %>%
  left_join(uniprot_gene, by = "gene") %>%
  left_join(go_gene, by = "gene") %>%
  left_join(surfaceome_gene, by = "gene")

gene_idx <- match(scored_tbl$gene, rownames(pct_mat))
state_idx <- match(scored_tbl$state, colnames(pct_mat))

scored_tbl$target_pct_cells <- pct_mat[cbind(gene_idx, state_idx)]
scored_tbl$target_mean_logexpr <- mean_mat[cbind(gene_idx, state_idx)]

get_off_target_summary <- function(mat, genes, states) {
  max_val <- rep(NA_real_, length(genes))
  max_state <- rep(NA_character_, length(genes))

  for (i in seq_along(genes)) {
    other_states <- setdiff(colnames(mat), states[i])
    vals <- mat[genes[i], other_states]
    if (length(vals) == 0 || all(is.na(vals))) {
      next
    }
    top_idx <- which.max(vals)
    max_val[i] <- vals[top_idx]
    max_state[i] <- other_states[top_idx]
  }

  list(value = max_val, state = max_state)
}

off_pct <- get_off_target_summary(pct_mat, scored_tbl$gene, scored_tbl$state)
off_expr <- get_off_target_summary(mean_mat, scored_tbl$gene, scored_tbl$state)

scored_tbl <- scored_tbl %>%
  mutate(
    max_off_target_pct_cells = off_pct$value,
    max_off_target_pct_state = off_pct$state,
    pct_margin = target_pct_cells - max_off_target_pct_cells,
    max_off_target_mean_logexpr = off_expr$value,
    max_off_target_expr_state = off_expr$state,
    expr_margin = target_mean_logexpr - max_off_target_mean_logexpr,
    go_surface_any_flag = coalesce(go_surface_flag, FALSE) | coalesce(go_plasma_membrane_flag, FALSE),
    annotation_surface_flag = coalesce(uniprot_surface_flag, FALSE) |
      go_surface_any_flag |
      coalesce(surfaceome_surface_flag, FALSE),
    extracellular_epitope_flag = coalesce(uniprot_extracellular_flag, FALSE) |
      coalesce(go_external_side_flag, FALSE) |
      coalesce(surfaceome_extracellular_flag, FALSE),
    membrane_anchor_flag = coalesce(uniprot_transmem_flag, FALSE) |
      coalesce(surfaceome_transmem_flag, FALSE),
    annotation_source_n = as.integer(coalesce(uniprot_surface_flag, FALSE)) +
      as.integer(go_surface_any_flag) +
      as.integer(coalesce(surfaceome_surface_flag, FALSE)),
    annotation_score = 2 * annotation_source_n +
      as.integer(extracellular_epitope_flag) +
      as.integer(membrane_anchor_flag) +
      as.integer(coalesce(cspa_support_flag, FALSE)),
    is_potential = annotation_surface_flag &
      extracellular_epitope_flag &
      (membrane_anchor_flag | coalesce(surfaceome_surface_flag, FALSE)) &
      pct_margin > 0 &
      expr_margin > 0 &
      (
        coalesce(surfaceome_surface_flag, FALSE) |
          (coalesce(uniprot_surface_flag, FALSE) & membrane_anchor_flag)
      ),
    is_supported = annotation_surface_flag &
      extracellular_epitope_flag &
      (membrane_anchor_flag | coalesce(surfaceome_surface_flag, FALSE)) &
      pct_margin > 0.1 &
      expr_margin > 0.5,
    supported_for_facs = is_potential | is_supported,
    recommended_for_facs = is_potential,
    facs_tier = case_when(
      is_potential ~ "potential",
      is_supported ~ "supported",
      TRUE ~ "weak_or_non_surface"
    )
  )

evidence_df <- scored_tbl %>%
  transmute(
    uni = coalesce(uniprot_surface_flag, FALSE),
    go = go_surface_any_flag,
    surfaceome = coalesce(surfaceome_surface_flag, FALSE),
    extracellular = extracellular_epitope_flag,
    transmem = membrane_anchor_flag,
    cspa = coalesce(cspa_support_flag, FALSE)
  )

scored_tbl$evidence_summary <- make_flag_summary(evidence_df)

scored_tbl <- scored_tbl %>%
  group_by(state) %>%
  mutate(
    annotation_rank = percent_rank0(annotation_score),
    target_pct_rank = percent_rank0(target_pct_cells),
    pct_margin_rank = percent_rank0(pct_margin),
    expr_margin_rank = percent_rank0(expr_margin),
    recurrence_rank = percent_rank0(sample_recurrence),
    facs_priority_score = as.integer(annotation_surface_flag) +
      2 * as.integer(recommended_for_facs) +
      3 * expr_margin_rank +
      3 * pct_margin_rank +
      target_pct_rank +
      recurrence_rank
  ) %>%
  ungroup() %>%
  arrange(
    state,
    desc(facs_priority_score),
    desc(recommended_for_facs),
    desc(supported_for_facs),
    desc(expr_margin),
    desc(pct_margin),
    desc(target_pct_cells),
    desc(sample_recurrence),
    gene
  ) %>%
  group_by(state) %>%
  mutate(state_rank = row_number()) %>%
  ungroup()

fwrite(
  scored_tbl,
  file.path(out_dir, "Auto_five_state_surface_marker_ranked.csv")
)

####################
# Excel workbook
####################
message("Writing ranked Excel workbook.")

summary_cols <- c(
  "state_rank",
  "state",
  "gene",
  "facs_tier",
  "target_pct_cells",
  "max_off_target_pct_cells",
  "target_mean_logexpr",
  "max_off_target_mean_logexpr",
  "sample_recurrence",
  "median_log2FC_hit",
  "evidence_summary"
)

all_scored_cols <- c(
  summary_cols,
  "uniprot_surface_flag",
  "go_surface_any_flag",
  "surfaceome_surface_flag",
  "extracellular_epitope_flag",
  "membrane_anchor_flag",
  "cspa_support_flag"
)

surface_export <- scored_tbl %>%
  filter(supported_for_facs) %>%
  dplyr::select(any_of(summary_cols))

all_export <- scored_tbl %>%
  dplyr::select(any_of(all_scored_cols))

wb <- createWorkbook()

top5_tbl <- scored_tbl %>%
  filter(supported_for_facs) %>%
  group_by(state) %>%
  slice_head(n = 5) %>%
  ungroup()

top5_export <- top5_tbl %>% dplyr::select(any_of(summary_cols))
top5_export$Sep_Main <- ""

expr_cols <- character()
for (st in state_order) {
  col_nm <- paste0("expr_", st)
  top5_export[[col_nm]] <- mean_mat[top5_export$gene, st]
  expr_cols <- c(expr_cols, col_nm)
  # Add separator
  sep_nm <- paste0("Sep_expr_", st)
  top5_export[[sep_nm]] <- ""
}

top5_export$Sep_Mid <- ""

pct_cols <- character()
for (st in state_order) {
  col_nm <- paste0("pct_", st)
  top5_export[[col_nm]] <- pct_mat[top5_export$gene, st]
  pct_cols <- c(pct_cols, col_nm)
  # Add separator
  sep_nm <- paste0("Sep_pct_", st)
  top5_export[[sep_nm]] <- ""
}

addWorksheet(wb, "Top5_per_state")
writeDataTable(wb, "Top5_per_state", top5_export)
freezePane(wb, "Top5_per_state", firstRow = TRUE)

addWorksheet(wb, "FACS_Candidates")
writeDataTable(wb, "FACS_Candidates", surface_export)
freezePane(wb, "FACS_Candidates", firstRow = TRUE)

addWorksheet(wb, "All_Scored")
writeDataTable(wb, "All_Scored", all_export)
freezePane(wb, "All_Scored", firstRow = TRUE)

for (state_name in state_order) {
  sheet_name <- safe_sheet_name(state_name)
  state_df <- scored_tbl %>%
    filter(state == state_name, supported_for_facs) %>%
    dplyr::select(any_of(summary_cols))

  addWorksheet(wb, sheet_name, tabColour = state_cols[state_name])
  writeDataTable(wb, sheet_name, state_df)
  freezePane(wb, sheet_name, firstRow = TRUE)
}

addWorksheet(wb, "DB_Manifest")
writeDataTable(wb, "DB_Manifest", db_manifest)
freezePane(wb, "DB_Manifest", firstRow = TRUE)

percent_style <- createStyle(numFmt = "0.0%")
score_style <- createStyle(numFmt = "0.000")
sep_style <- createStyle(fgFill = "#D5D8DC", border = "LeftRight", borderColour = "#95A5A6")

percent_cols <- c(
  "target_pct_cells",
  "max_off_target_pct_cells",
  "sample_recurrence",
  pct_cols
)

score_cols <- c(
  "target_mean_logexpr",
  "max_off_target_mean_logexpr",
  "median_log2FC_hit",
  expr_cols
)

format_sheet <- function(sheet_name, df) {
  if (nrow(df) == 0) return(invisible(NULL))

  # Apply state-specific coloring to the 'state' column
  state_col_idx <- match("state", names(df))
  if (!is.na(state_col_idx)) {
    for (st_name in names(state_cols)) {
      st_style <- createStyle(fontColour = "#FFFFFF", fgFill = state_cols[st_name], textDecoration = "bold")
      rows <- which(df$state == st_name) + 1
      if (length(rows) > 0) {
        addStyle(wb, sheet = sheet_name, style = st_style, rows = rows, cols = state_col_idx, gridExpand = TRUE, stack = TRUE)
      }
    }
  }

  # Apply state-specific coloring to expression and pct headers in Top5 sheet
  if (sheet_name == "Top5_per_state") {
    for (st_name in names(state_cols)) {
      h_style <- createStyle(fontColour = "#FFFFFF", fgFill = state_cols[st_name], textDecoration = "bold", halign = "center")
      
      e_col <- match(paste0("expr_", st_name), names(df))
      if (!is.na(e_col)) {
        addStyle(wb, sheet = sheet_name, style = h_style, rows = 1, cols = e_col, stack = TRUE)
      }
      
      p_col <- match(paste0("pct_", st_name), names(df))
      if (!is.na(p_col)) {
        addStyle(wb, sheet = sheet_name, style = h_style, rows = 1, cols = p_col, stack = TRUE)
      }
    }
  }

  for (col_name in percent_cols) {
    idx <- match(col_name, names(df))
    if (!is.na(idx)) {
      addStyle(
        wb,
        sheet = sheet_name,
        style = percent_style,
        rows = 2:(nrow(df) + 1),
        cols = idx,
        gridExpand = TRUE,
        stack = TRUE
      )
      
      # Colour formatting for pct
      conditionalFormatting(
        wb, sheet_name, cols = idx, rows = 2:(nrow(df) + 1),
        style = c("#FFFFFF", "#27AE60"),
        rule = c(0, 1),
        type = "colourScale"
      )
    }
  }

  for (col_name in score_cols) {
    idx <- match(col_name, names(df))
    if (!is.na(idx)) {
      addStyle(
        wb,
        sheet = sheet_name,
        style = score_style,
        rows = 2:(nrow(df) + 1),
        cols = idx,
        gridExpand = TRUE,
        stack = TRUE
      )
      
      # Colour formatting for expr
      conditionalFormatting(
        wb, sheet_name, cols = idx, rows = 2:(nrow(df) + 1),
        style = c("#FFFFFF", "#FB8A8A", "#B22222"),
        rule = c(0, 2.5, 5),
        type = "colourScale"
      )
    }
  }

  # Add sep style
  sep_cols <- names(df)[grepl("^Sep", names(df))]
  for (col_name in sep_cols) {
    idx <- match(col_name, names(df))
    if (!is.na(idx)) {
      addStyle(wb, sheet = sheet_name, style = sep_style, rows = 1:(nrow(df) + 1), cols = idx, stack = TRUE)
      setColWidths(wb, sheet = sheet_name, cols = idx, widths = 2)
    }
  }

  non_sep_cols <- which(!names(df) %in% sep_cols)
  setColWidths(wb, sheet = sheet_name, cols = non_sep_cols, widths = "auto")
}

format_sheet("Top5_per_state", top5_export)
format_sheet("FACS_Candidates", surface_export)
format_sheet("All_Scored", all_export)

for (state_name in state_order) {
  sheet_name <- safe_sheet_name(state_name)
  state_df <- scored_tbl %>%
    filter(state == state_name, supported_for_facs) %>%
    dplyr::select(any_of(summary_cols))
  format_sheet(sheet_name, state_df)
}

setColWidths(wb, sheet = "DB_Manifest", cols = 1:ncol(db_manifest), widths = "auto")

saveWorkbook(
  wb,
  file = file.path(out_dir, "Auto_five_state_surface_marker_candidates.xlsx"),
  overwrite = TRUE
)

message("Generating bubble plot for top 5 markers...")
library(ggplot2)

plot_data <- data.frame(
  gene = rep(top5_export$gene, length(state_order)),
  target_state = rep(top5_export$state, length(state_order)),
  state = rep(state_order, each = nrow(top5_export)),
  pct = as.vector(as.matrix(pct_data)),
  expr = as.vector(as.matrix(expr_data))
)
plot_data$gene <- factor(plot_data$gene, levels = rev(unique(top5_export$gene)))
plot_data$state <- factor(plot_data$state, levels = state_order)
plot_data$target_state <- factor(plot_data$target_state, levels = state_order)

p <- ggplot(plot_data, aes(x = state, y = gene, size = pct, color = expr)) +
  geom_point() +
  scale_size_continuous(range = c(0, 6)) +
  scale_color_gradient(low = "lightgrey", high = "red", limits = c(0, max(plot_data$expr, na.rm=TRUE))) +
  facet_grid(target_state ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0),
    panel.grid.major = element_line(color = "gray90")
  ) +
  labs(title = "Top 5 Surface Markers Expression", x = "State", y = "Gene", size = "Percent Expressed", color = "Median Mean Expr")

ggsave(file.path(out_dir, "Auto_five_state_surface_marker_dotplot.pdf"), plot = p, width = 8, height = max(6, length(unique(top5_export$gene)) * 0.3), useDingbats = FALSE)

message("Finished surface-marker prioritization.")
