####################
# Auto_drug_reversal_scdrugprio.R
#
# Run scDrugPrio transcriptomic-reversal drug prioritization from PDO DEGs.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

####################
# setup
####################

project_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
setwd(file.path(project_dir, "PDOs_outs"))

base_dir <- "Auto_drug_reversal"
input_dir <- file.path(base_dir, "scdrugprio_inputs")
out_dir <- file.path(base_dir, "scdrugprio")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

state_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation"
)

params <- list(
  cores = as.integer(Sys.getenv("AUTO_SCDRUGPRIO_CORES", "8")),
  n_random_iterations = as.integer(Sys.getenv("AUTO_SCDRUGPRIO_RANDOM_ITERATIONS", "1000")),
  p_adj_threshold = as.numeric(Sys.getenv("AUTO_SCDRUGPRIO_PADJ", "0.05")),
  min_abs_logfc = as.numeric(Sys.getenv("AUTO_SCDRUGPRIO_MIN_ABS_LOGFC", "0.10")),
  max_disease_genes = as.integer(Sys.getenv("AUTO_SCDRUGPRIO_MAX_DISEASE_GENES", "1500")),
  exclude_mimics = !identical(Sys.getenv("AUTO_SCDRUGPRIO_EXCLUDE_MIMICS", "1"), "0"),
  require_counteracting_targets = !identical(Sys.getenv("AUTO_SCDRUGPRIO_REQUIRE_COUNTERACTION", "1"), "0"),
  reuse_state_results = !identical(Sys.getenv("AUTO_SCDRUGPRIO_REUSE_STATE_RESULTS", "1"), "0")
)

####################
# helpers
####################

safe_state_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

write_status <- function(status, detail) {
  fwrite(
    data.frame(
      step = "scdrugprio",
      status = status,
      detail = detail,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Auto_scdrugprio_status.csv")
  )
}

first_matching_col <- function(df, candidates) {
  nms <- colnames(df)
  hit <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit) > 0) return(hit[1])
  pattern <- paste(candidates, collapse = "|")
  hit <- grep(pattern, nms, ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) hit[1] else NA_character_
}

read_matrix_file <- function(path) {
  x <- fread(path)
  as.matrix(x)
}

read_resource_table <- function(path, required = 2) {
  x <- fread(path)
  if (ncol(x) < required) stop("Resource file has too few columns: ", path)
  x
}

load_default_resource <- function(resource_name) {
  local_path <- file.path(base_dir, "resources", paste0("Auto_", resource_name, ".rda"))
  if (file.exists(local_path)) {
    load(local_path, envir = environment())
    return(get(resource_name, envir = environment()))
  }
  data(list = resource_name, package = "scDrugPrio", envir = environment())
  get(resource_name, envir = environment())
}

prepare_scdrugprio_functions <- function() {
  source_dir <- file.path(base_dir, "resources", "scdrugprio_source")
  dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

  files <- c(
    "average_closest_distance.R",
    "bin_creation_by_min_bin_size.R",
    "random_drug_target_bin_adjusted_distances.R",
    "create_drug_target_matrix.R",
    "average_closest_distance_network_drug_screening.R"
  )

  base_url <- "https://raw.githubusercontent.com/SDTC-CPMed/scDrugPrio/main/R"
  for (f in files) {
    dest <- file.path(source_dir, f)
    if (!file.exists(dest)) {
      download.file(paste(base_url, f, sep = "/"), destfile = dest, quiet = TRUE)
    }
  }

  patch_source <- function(path) {
    txt <- readLines(path, warn = FALSE)
    txt <- gsub('.packages = c("scDrugPrio", "doParallel")', '.packages = c("doParallel", "igraph")', txt, fixed = TRUE)
    txt <- gsub('.packages = "scDrugPrio"', '.packages = character()', txt, fixed = TRUE)
    eval(parse(text = txt), envir = .GlobalEnv)
  }

  invisible(lapply(file.path(source_dir, files), patch_source))

  list(
    create_drug_target_matrix = get("create_drug_target_matrix", envir = .GlobalEnv),
    average_closest_distance_network_drug_screening = get("average_closest_distance_network_drug_screening", envir = .GlobalEnv)
  )
}

map_symbol_targets_to_ppi_space <- function(ppi, drug_targets, target_col, deg_tables) {
  ppi_values <- unique(as.character(c(ppi[, 1], ppi[, 2])))
  ppi_values <- ppi_values[!is.na(ppi_values) & nzchar(ppi_values)]
  numeric_ppi <- mean(grepl("^[0-9]+$", ppi_values[seq_len(min(length(ppi_values), 5000))])) > 0.90

  if (!numeric_ppi) {
    return(list(
      drug_targets = drug_targets,
      deg_tables = deg_tables,
      id_space = "gene_symbol",
      symbol_map = NULL
    ))
  }

  if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
      !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("PPI appears to use Entrez IDs; AnnotationDbi/org.Hs.eg.db is required for symbol mapping.")
  }

  all_symbols <- unique(c(
    as.character(drug_targets[[target_col]]),
    unlist(lapply(deg_tables, function(x) x$gene))
  ))
  all_symbols <- all_symbols[!is.na(all_symbols) & nzchar(all_symbols)]

  mapped <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = all_symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  mapped <- mapped[!is.na(mapped$ENTREZID) & !duplicated(mapped$SYMBOL), , drop = FALSE]
  symbol_map <- mapped$ENTREZID
  names(symbol_map) <- mapped$SYMBOL

  drug_targets[[target_col]] <- unname(symbol_map[as.character(drug_targets[[target_col]])])
  drug_targets <- drug_targets[!is.na(drug_targets[[target_col]]) & nzchar(drug_targets[[target_col]]), , drop = FALSE]

  deg_tables <- lapply(deg_tables, function(df) {
    df$gene_symbol <- df$gene
    df$gene <- unname(symbol_map[as.character(df$gene_symbol)])
    df <- df[!is.na(df$gene) & nzchar(df$gene), , drop = FALSE]
    df
  })

  list(
    drug_targets = drug_targets,
    deg_tables = deg_tables,
    id_space = "entrez_id",
    symbol_map = symbol_map
  )
}

infer_action_col <- function(df) {
  nms <- colnames(df)
  action_candidates <- c("drug_action", "action", "pharmacological_action", "effect", "target_organism")
  hits <- nms[tolower(nms) %in% tolower(action_candidates)]
  action_words <- "agon|antagon|inhib|activ|block|suppress|induc|stimulat|potentiat|binder|cofactor|modulat"

  scored <- vapply(nms, function(nm) {
    values <- tolower(as.character(df[[nm]]))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) == 0) return(0)
    mean(grepl(action_words, values[seq_len(min(length(values), 5000))]))
  }, numeric(1))

  if (length(hits) > 0 && max(scored[hits], na.rm = TRUE) > 0.05) {
    return(hits[which.max(scored[hits])])
  }
  if (max(scored, na.rm = TRUE) > 0.05) names(scored)[which.max(scored)] else NA_character_
}

standardize_pharma_effect <- function(path, fallback_targets, fallback_drug_col, fallback_target_col, symbol_map = NULL) {
  if (file.exists(path)) {
    pharma_raw <- read_resource_table(path, required = 3)
  } else {
    pharma_raw <- fallback_targets
  }

  drug_col <- first_matching_col(pharma_raw, c("drug", "drug_id", "DrugBankID", "drugID", "pert_iname", "compound"))
  target_col <- first_matching_col(pharma_raw, c("gene_symbol", "target_gene", "target", "gene", "Gene", "name"))
  action_col <- infer_action_col(pharma_raw)

  if (is.na(drug_col)) drug_col <- fallback_drug_col
  if (is.na(target_col)) target_col <- fallback_target_col
  if (is.na(action_col) || is.na(drug_col) || is.na(target_col)) return(NULL)

  out <- tibble(
    drug = as.character(pharma_raw[[drug_col]]),
    target = as.character(pharma_raw[[target_col]]),
    action = as.character(pharma_raw[[action_col]])
  ) %>%
    filter(!is.na(drug), nzchar(drug), !is.na(target), nzchar(target), !is.na(action), nzchar(action))

  if (!is.null(symbol_map)) {
    out$target <- unname(symbol_map[out$target])
    out <- out %>% filter(!is.na(target), nzchar(target))
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}

counteraction_summary <- function(drug_ids, deg_df, pharma_effect) {
  if (is.null(pharma_effect) || nrow(pharma_effect) == 0) {
    return(tibble(
      drug = drug_ids,
      counteracting_targets = NA_character_,
      mimicking_targets = NA_character_,
        n_counteracting_targets = NA_integer_,
        n_mimicking_targets = NA_integer_,
        n_directional_targets = NA_integer_,
        counteract_fraction = NA_real_,
        mimic_fraction = NA_real_,
        net_direction_score = NA_real_,
        direction_call = "no_action_resource"
      ))
  }

  drug_col <- colnames(pharma_effect)[1]
  target_col <- colnames(pharma_effect)[2]
  action_col <- colnames(pharma_effect)[3]
  deg_map <- deg_df$avg_logFC
  names(deg_map) <- deg_df$gene

  bind_rows(lapply(drug_ids, function(drug_id) {
    rows <- pharma_effect[as.character(pharma_effect[[drug_col]]) == drug_id, , drop = FALSE]
    if (nrow(rows) == 0) {
      return(tibble(
        drug = drug_id,
        counteracting_targets = NA_character_,
        mimicking_targets = NA_character_,
        n_counteracting_targets = 0L,
        n_mimicking_targets = 0L,
        n_directional_targets = 0L,
        counteract_fraction = 0,
        mimic_fraction = 0,
        net_direction_score = 0,
        direction_call = "no_directional_target"
      ))
    }

    targets <- unique(as.character(rows[[target_col]]))
    actions <- rows[[action_col]][match(targets, as.character(rows[[target_col]]))]
    logfc <- deg_map[targets]
    action_lower <- tolower(as.character(actions))
    activates <- grepl("agon|activ|stimulat|induc|posit", action_lower)
    inhibits <- grepl("antagon|inhib|block|suppress|negative|inverse", action_lower)

    counteracts <- (!is.na(logfc) & logfc > 0 & inhibits) | (!is.na(logfc) & logfc < 0 & activates)
    mimics <- (!is.na(logfc) & logfc > 0 & activates) | (!is.na(logfc) & logfc < 0 & inhibits)

    counteracting_targets <- targets[which(counteracts)]
    mimicking_targets <- targets[which(mimics)]
    n_targeted_degs <- length(which(!is.na(logfc)))
    n_mimicking <- length(mimicking_targets)
    n_counteracting <- length(counteracting_targets)
    direction_call <- dplyr::case_when(
      n_targeted_degs == 0 ~ "no_directional_target",
      n_counteracting > 0 & n_mimicking == 0 ~ "counteracting",
      n_counteracting > 0 & n_mimicking > 0 ~ "mixed",
      n_counteracting == 0 & n_mimicking > 0 ~ "mimicking",
      TRUE ~ "non_directional_action"
    )

    tibble(
      drug = drug_id,
      counteracting_targets = paste(counteracting_targets, collapse = ";"),
      mimicking_targets = paste(mimicking_targets, collapse = ";"),
      n_counteracting_targets = n_counteracting,
      n_mimicking_targets = n_mimicking,
      n_directional_targets = n_targeted_degs,
      counteract_fraction = ifelse(n_targeted_degs > 0, n_counteracting / n_targeted_degs, 0),
      mimic_fraction = ifelse(n_targeted_degs > 0, n_mimicking / n_targeted_degs, 0),
      net_direction_score = ifelse(n_targeted_degs > 0, (n_counteracting - n_mimicking) / n_targeted_degs, 0),
      direction_call = direction_call
    )
  }))
}

normalize_scdrug_table <- function(df) {
  df <- as_tibble(df)
  char_cols <- intersect(
    c(
      "state", "pipeline", "drug", "drug_id", "direction_call", "target_genes",
      "counteracting_targets", "mimicking_targets"
    ),
    colnames(df)
  )
  num_cols <- intersect(
    c(
      "rank", "network_rank", "score", "p_value", "dc", "zc", "n_drug_targets",
      "direction_priority", "net_direction_score", "counteract_fraction",
      "mimic_fraction", "n_directional_targets", "n_counteracting_targets",
      "n_mimicking_targets"
    ),
    colnames(df)
  )
  df <- df %>% mutate(across(all_of(char_cols), as.character))
  df <- df %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))
  df
}

####################
# validation
####################

use_package <- requireNamespace("scDrugPrio", quietly = TRUE)
scdrug_fns <- NULL
if (!use_package) {
  scdrug_fns <- tryCatch(
    prepare_scdrugprio_functions(),
    error = function(e) {
      write_status("missing_package", paste("scDrugPrio package is absent and source fallback failed:", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(scdrug_fns)) quit(save = "no", status = 0)
}

ppi_path <- Sys.getenv("AUTO_SCDRUGPRIO_PPI", "")
drug_targets_path <- Sys.getenv("AUTO_SCDRUGPRIO_DRUG_TARGETS", "")
pharma_effect_path <- Sys.getenv("AUTO_SCDRUGPRIO_PHARMA_EFFECT", "")

deg_files <- file.path(input_dir, paste0("Auto_scdrugprio_deg_", safe_state_name(state_order), ".txt"))
if (any(!file.exists(deg_files))) {
  write_status("missing_inputs", "Missing scDrugPrio DEG files. Run Auto_drug_reversal_inputs.R first.")
  quit(save = "no", status = 0)
}

####################
# resources
####################

message("Loading scDrugPrio resources.")
used_default_resources <- FALSE
if (file.exists(ppi_path)) {
  ppi <- read_matrix_file(ppi_path)[, 1:2, drop = FALSE]
} else {
  ppi <- load_default_resource("lit_ppi")
  used_default_resources <- TRUE
}
ppi <- matrix(as.character(ppi), ncol = 2, dimnames = dimnames(ppi))

if (file.exists(drug_targets_path)) {
  drug_targets <- read_resource_table(drug_targets_path, required = 2)
} else {
  drug_targets <- as.data.frame(load_default_resource("drug_bank_example_data"), stringsAsFactors = FALSE)
  used_default_resources <- TRUE
}
raw_drug_targets <- drug_targets

drug_col <- first_matching_col(drug_targets, c("drug", "drug_id", "DrugBankID", "drugID", "pert_iname", "compound"))
target_col <- first_matching_col(drug_targets, c("gene_symbol", "target_gene", "target", "gene", "Gene", "name"))
drug_name_col <- first_matching_col(drug_targets, c("drug_name", "pert_iname", "compound_name", "name"))
if (is.na(drug_col)) drug_col <- colnames(drug_targets)[1]
if (is.na(target_col)) target_col <- colnames(drug_targets)[2]
drug_display_map <- setNames(as.character(drug_targets[[if (!is.na(drug_name_col)) drug_name_col else drug_col]]), as.character(drug_targets[[drug_col]]))
drug_display_map <- drug_display_map[!duplicated(names(drug_display_map))]

deg_tables <- lapply(deg_files, fread)
names(deg_tables) <- state_order

mapped_inputs <- map_symbol_targets_to_ppi_space(
  ppi = ppi,
  drug_targets = drug_targets,
  target_col = target_col,
  deg_tables = deg_tables
)
drug_targets <- mapped_inputs$drug_targets
deg_tables <- mapped_inputs$deg_tables

create_drug_target_matrix_fn <- if (use_package) scDrugPrio::create_drug_target_matrix else scdrug_fns$create_drug_target_matrix
average_distance_fn <- if (use_package) scDrugPrio::average_closest_distance_network_drug_screening else scdrug_fns$average_closest_distance_network_drug_screening

drug_target_matrix <- create_drug_target_matrix_fn(
  drugID = as.character(drug_targets[[drug_col]]),
  target = as.character(drug_targets[[target_col]])
)

pharma_effect <- standardize_pharma_effect(
  path = pharma_effect_path,
  fallback_targets = raw_drug_targets,
  fallback_drug_col = drug_col,
  fallback_target_col = target_col,
  symbol_map = mapped_inputs$symbol_map
)

ranked_all <- list()

####################
# state-wise prioritization
####################

for (state_name in state_order) {
  state_safe <- safe_state_name(state_name)
  state_dir <- file.path(out_dir, state_safe)
  dir.create(state_dir, recursive = TRUE, showWarnings = FALSE)
  ranked_path <- file.path(state_dir, paste0("Auto_scdrugprio_ranked_", state_safe, ".csv"))
  audit_path <- file.path(state_dir, paste0("Auto_scdrugprio_direction_audit_", state_safe, ".csv"))

  if (params$reuse_state_results && file.exists(ranked_path) && file.exists(audit_path)) {
    message("Reusing completed scDrugPrio state result: ", state_name)
    ranked_all[[state_name]] <- normalize_scdrug_table(fread(ranked_path))
    next
  }

  deg_df <- deg_tables[[state_name]]
  disease_genes <- deg_df %>%
    filter(p_val_adj <= params$p_adj_threshold, abs(avg_logFC) >= params$min_abs_logfc) %>%
    arrange(p_val_adj, desc(abs(avg_logFC))) %>%
    pull(gene) %>%
    unique()

  if (length(disease_genes) < 10) {
    disease_genes <- deg_df %>%
      arrange(p_val_adj, desc(abs(avg_logFC))) %>%
      slice_head(n = 300) %>%
      pull(gene) %>%
      unique()
  }

  if (is.finite(params$max_disease_genes) && params$max_disease_genes > 0 && length(disease_genes) > params$max_disease_genes) {
    disease_genes <- disease_genes[seq_len(params$max_disease_genes)]
  }

  message("Running scDrugPrio network screen: ", state_name, " (", length(disease_genes), " disease genes)")

  dist_res <- average_distance_fn(
    ppin = ppi,
    drug_target_matrix = drug_target_matrix,
    disease_genes = disease_genes,
    file_name = paste0("Auto_", state_safe),
    disease_genes_lcc = FALSE,
    n_random_iterations = params$n_random_iterations,
    cores = params$cores,
    out_dir = state_dir
  )

  dist_df <- as.data.frame(dist_res, stringsAsFactors = FALSE)
  colnames(dist_df) <- make.names(colnames(dist_df))
  dist_df <- dist_df %>%
    transmute(
      drug_id = as.character(Drug),
      drug = unname(drug_display_map[drug_id]),
      drug = ifelse(is.na(drug) | !nzchar(drug), drug_id, drug),
      n_drug_targets = suppressWarnings(as.numeric(n_drug_targets)),
      dc = suppressWarnings(as.numeric(dc)),
      zc = suppressWarnings(as.numeric(zc)),
      p_value = suppressWarnings(as.numeric(P))
    )

  selected <- dist_df %>%
    filter(is.finite(p_value), p_value < 0.05, is.finite(dc), dc < 1)

  if (nrow(selected) == 0) {
    selected <- dist_df %>%
      filter(is.finite(p_value), is.finite(dc)) %>%
      arrange(p_value, dc) %>%
      slice_head(n = 500)
  }

  counter <- counteraction_summary(selected$drug_id, deg_df, pharma_effect)

  ranked_candidates <- selected %>%
    left_join(counter, by = c("drug_id" = "drug")) %>%
    mutate(
      state = state_name,
      pipeline = "scDrugPrio",
      direction_priority = case_when(
        direction_call == "counteracting" ~ 1L,
        direction_call == "mixed" ~ 2L,
        direction_call %in% c("no_directional_target", "non_directional_action", "no_action_resource") ~ 3L,
        direction_call == "mimicking" ~ 4L,
        TRUE ~ 5L
      ),
      score = -zc + ifelse(is.na(net_direction_score), 0, net_direction_score),
      target_genes = ifelse(
        !is.na(counteracting_targets) & nzchar(counteracting_targets),
        counteracting_targets,
        NA_character_
      )
    )

  ranked_audit <- ranked_candidates %>%
    arrange(p_value, dc, desc(net_direction_score), desc(counteract_fraction), drug) %>%
    mutate(network_rank = row_number()) %>%
    select(
      state, pipeline, drug, drug_id, network_rank, score, p_value, dc, zc, n_drug_targets,
      direction_call, direction_priority, net_direction_score,
      counteract_fraction, mimic_fraction, n_directional_targets,
      n_counteracting_targets, n_mimicking_targets,
      target_genes, counteracting_targets, mimicking_targets
    )
  ranked_audit <- normalize_scdrug_table(ranked_audit)

  fwrite(ranked_audit, audit_path)

  if (params$exclude_mimics && any(ranked_candidates$direction_call != "mimicking", na.rm = TRUE)) {
    ranked_candidates <- ranked_candidates %>% filter(direction_call != "mimicking")
  }

  if (params$require_counteracting_targets) {
    ranked_candidates <- ranked_candidates %>% filter(n_counteracting_targets > 0, net_direction_score > 0)
  }

  ranked <- ranked_candidates %>%
    arrange(
      direction_priority,
      desc(net_direction_score),
      desc(counteract_fraction),
      p_value,
      dc,
      desc(score),
      drug
    ) %>%
    mutate(rank = row_number()) %>%
    select(
      state, pipeline, drug, drug_id, rank, score, p_value, dc, zc, n_drug_targets,
      direction_call, direction_priority, net_direction_score,
      counteract_fraction, mimic_fraction, n_directional_targets,
      n_counteracting_targets, n_mimicking_targets,
      target_genes, counteracting_targets, mimicking_targets
    )
  ranked <- normalize_scdrug_table(ranked)

  fwrite(ranked, ranked_path)
  ranked_all[[state_name]] <- ranked
}

ranked_all <- bind_rows(lapply(ranked_all, normalize_scdrug_table))
audit_all <- bind_rows(lapply(state_order, function(state_name) {
  state_safe <- safe_state_name(state_name)
  path <- file.path(out_dir, state_safe, paste0("Auto_scdrugprio_direction_audit_", state_safe, ".csv"))
  if (file.exists(path)) normalize_scdrug_table(fread(path)) else NULL
}))
if (nrow(audit_all) > 0) {
  fwrite(audit_all, file.path(out_dir, "Auto_scdrugprio_direction_audit_all_drugs.csv"))
}

if (nrow(ranked_all) == 0) {
  write_status("no_results", "scDrugPrio completed but no rankable drug rows were returned.")
  quit(save = "no", status = 0)
}

fwrite(ranked_all, file.path(out_dir, "Auto_scdrugprio_ranked_drugs.csv"))
write_status(
  "complete",
  paste(
    "Ranked", nrow(ranked_all), "scDrugPrio drug-state rows.",
    "Resource mode:", ifelse(used_default_resources, "package_default_lit_ppi_drug_bank_example_data", "explicit_user_resources"),
    "Function source:", ifelse(use_package, "installed_scDrugPrio_package", "sourced_scDrugPrio_R_functions"),
    "Identifier space:", mapped_inputs$id_space,
    "Max disease genes:", params$max_disease_genes,
    "Exclude mimicking target actions:", params$exclude_mimics,
    "Require counteracting target evidence:", params$require_counteracting_targets
  )
)

message("scDrugPrio drug prioritization complete.")
