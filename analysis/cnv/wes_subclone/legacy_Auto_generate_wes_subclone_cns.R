####################
# Analysis registry (authoritative override):
#   Status: legacy; retained for provenance, no current downstream use
#   Script: analysis/cnv/wes_subclone/legacy_Auto_generate_wes_subclone_cns.R
#   Methodology: historical method only; see analysis/ANALYSIS_MAP.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description:
#     Preserves a superseded implementation or analysis tied to superseded
#     inputs. Do not use its outputs as current centred-MP/state inputs. The
#     original historical inputs, outputs, and method notes remain below.
####################

####################
# Auto_generate_wes_subclone_cns.R
#
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_generate_wes_subclone_cns.R
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_segments.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_purity_ploidy.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_results.tsv
#   - PDOs_outs/Auto_wes_subclone/tables/pyclone/Auto_<sample>_pyclone_vi_mutation_map.tsv
#   - sarek_mutect CNVkit .cns (whole-sample bulk reference)
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/tables/cns/Auto_<sample>_bulk.cns
#   - PDOs_outs/Auto_wes_subclone/tables/cns/Auto_<sample>_cluster<N>.cns
#   - PDOs_outs/Auto_wes_subclone/tables/cns/Auto_wes_subclone_cns_summary.csv
# Downstream use: per-subclone CNA profiles for cross-platform clone matching.
####################

suppressPackageStartupMessages({
  library(data.table)
})

# ── Configuration ─────────────────────────────────────────────────────────────
out_root <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone"
cnvkit_root <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/cnvkit"

samples <- c("PDO_1090_vs_NT_1090", "PDO_1181_vs_NT_1181")

cns_dir <- file.path(out_root, "tables/cns")
dir.create(cns_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helper: compute log2 ratio from total copy number and ploidy ──────────────
# For a subclone at given CCF carrying a segment at total_cn, the expected
# log2 ratio relative to diploid (adjusted for purity) is:
#   log2( (purity * total_cn + (1-purity) * 2) / (purity * ploidy + (1-purity) * 2) )
compute_log2 <- function(total_cn, purity, ploidy) {
  expected <- purity * total_cn + (1 - purity) * 2
  baseline <- purity * ploidy + (1 - purity) * 2
  log2(expected / baseline)
}

# ── Main loop per sample ─────────────────────────────────────────────────────
summary_rows <- list()

for (sample in samples) {
  message("\n=== Processing ", sample, " ===")

  # ── Load FACETS segments ──────────────────────────────────────────────────
  seg_path <- file.path(out_root, "tables/facets",
                        paste0("Auto_", sample, "_facets_segments.tsv"))
  pp_path <- file.path(out_root, "tables/facets",
                       paste0("Auto_", sample, "_facets_purity_ploidy.tsv"))
  if (!file.exists(seg_path) || !file.exists(pp_path)) {
    message("  Missing FACETS tables, skipping ", sample)
    next
  }
  segments <- fread(seg_path)
  pp <- fread(pp_path)
  purity <- as.numeric(pp$purity[1])
  ploidy <- as.numeric(pp$ploidy[1])
  message("  FACETS purity=", round(purity, 3), " ploidy=", round(ploidy, 3))

  # ── Load PyClone-VI results + mutation map ────────────────────────────────
  res_path <- file.path(out_root, "tables/pyclone",
                        paste0("Auto_", sample, "_pyclone_vi_results.tsv"))
  map_path <- file.path(out_root, "tables/pyclone",
                        paste0("Auto_", sample, "_pyclone_vi_mutation_map.tsv"))
  if (!file.exists(res_path) || !file.exists(map_path)) {
    message("  Missing PyClone tables, skipping ", sample)
    next
  }
  pyclone_res <- fread(res_path)
  pyclone_map <- fread(map_path)

  # Merge cluster assignments onto mutation positions
  pyclone_res[, cluster_id := as.character(cluster_id)]
  merged <- merge(pyclone_res[, .(mutation_id, cluster_id, cellular_prevalence,
                                   cluster_assignment_prob)],
                  pyclone_map[, .(mutation_id, chrom, pos)],
                  by = "mutation_id")
  merged[, chrom_clean := sub("^chr", "", as.character(chrom))]
  merged[, pos := as.integer(pos)]

  cluster_ids <- sort(unique(merged$cluster_id))
  n_clusters <- length(cluster_ids)
  message("  PyClone clusters: ", paste(cluster_ids, collapse = ", "),
          " (", n_clusters, " total)")

  # ── Compute per-cluster median CCF ────────────────────────────────────────
  cluster_ccf <- merged[, .(median_ccf = median(cellular_prevalence, na.rm = TRUE),
                            n_mutations = .N),
                        by = cluster_id]
  setorder(cluster_ccf, -median_ccf)
  message("  Cluster CCFs: ",
          paste0("c", cluster_ccf$cluster_id, "=",
                 round(cluster_ccf$median_ccf, 3), " (n=", cluster_ccf$n_mutations, ")"),
          collapse = ", ")

  # ── Prepare segments for .cns output ──────────────────────────────────────
  # Standardise chromosome naming to match CNVkit format (chr1, chr2, ...)
  segments[, chrom_clean := sub("^chr", "", as.character(chrom))]
  autosomes <- as.character(1:22)
  segments <- segments[chrom_clean %in% autosomes]
  segments[, chromosome := paste0("chr", chrom_clean)]
  segments[, start := as.integer(start)]
  segments[, end := as.integer(end)]

  # Compute bulk log2 for each segment
  segments[, log2_bulk := compute_log2(total_cn, purity, ploidy)]

  # ── 1. Write bulk .cns ──────────────────────────────────────────────────
  bulk_cns <- segments[, .(
    chromosome, start, end,
    gene = "-",
    log2 = round(log2_bulk, 6),
    depth = 0,
    probes = 0,
    weight = 1
  )]

  # Load the original CNVkit .cns for a more accurate high-resolution bulk profile
  tumor_id <- sub("_vs_.*", "", sample)
  cnvkit_cns_path <- file.path(cnvkit_root, sample, paste0(tumor_id, ".cns"))
  if (file.exists(cnvkit_cns_path)) {
    message("  Using original CNVkit .cns as bulk reference")
    cnvkit_cns <- fread(cnvkit_cns_path)
    # Write as-is in standard format
    bulk_cns <- cnvkit_cns[, .(
      chromosome, start, end,
      gene = if ("gene" %in% names(cnvkit_cns)) cnvkit_cns$gene else "-",
      log2,
      depth = if ("depth" %in% names(cnvkit_cns)) cnvkit_cns$depth else 0,
      probes = if ("probes" %in% names(cnvkit_cns)) cnvkit_cns$probes else 0,
      weight = if ("weight" %in% names(cnvkit_cns)) cnvkit_cns$weight else 1
    )]
  }

  bulk_path <- file.path(cns_dir, paste0("Auto_", sample, "_bulk.cns"))
  fwrite(bulk_cns, bulk_path, sep = "\t")
  message("  Wrote bulk .cns: ", bulk_path, " (", nrow(bulk_cns), " segments)")

  # ── 2. Per-subclone .cns files ────────────────────────────────────────────
  # Strategy: For each cluster, we use the FACETS cellular fraction (cf) to 
  # determine if the subclone possesses the CNA.
  # If cf >= cluster_ccf - 0.15, the CNA is present in at least as many cells 
  # as the subclone, so the subclone inherits it.
  # Otherwise, the CNA is in a smaller subpopulation, so this subclone does not have it.

  for (cid in cluster_ids) {
    cluster_median_ccf <- cluster_ccf[cluster_id == cid, median_ccf]

    cluster_cns <- copy(segments)
    
    # Assign log2 based on cf (scaled to cancer cell fraction by dividing by purity)
    cluster_cns[, log2_subclone := fifelse(
      is.na(cf) | (cf / purity) >= (cluster_median_ccf - 0.15),
      log2_bulk,
      0
    )]

    # Write .cns
    out_cns <- cluster_cns[, .(
      chromosome, start, end,
      gene = "-",
      log2 = round(log2_subclone, 6),
      depth = 0,
      probes = fifelse(log2_subclone != 0, 1, 0),
      weight = 1
    )]

    cns_path <- file.path(cns_dir,
                          paste0("Auto_", sample, "_cluster", cid, ".cns"))
    fwrite(out_cns, cns_path, sep = "\t")
    message("  Wrote cluster ", cid, " .cns: ", cns_path,
            " (CCF=", round(cluster_median_ccf, 3), ")")

    summary_rows[[paste0(sample, "_", cid)]] <- data.table(
      sample = sample,
      cluster_id = cid,
      median_ccf = cluster_median_ccf,
      n_mutations = cluster_ccf[cluster_id == cid, n_mutations],
      n_segments = nrow(cluster_cns),
      n_segments_with_muts = sum(cluster_cns$log2_subclone != 0),
      cns_path = cns_path,
      purity = purity,
      ploidy = ploidy
    )
  }
}

# ── Summary table ─────────────────────────────────────────────────────────────
if (length(summary_rows) > 0) {
  summary_dt <- rbindlist(summary_rows)
  summary_path <- file.path(cns_dir, "Auto_wes_subclone_cns_summary.csv")
  fwrite(summary_dt, summary_path)
  message("\nWrote .cns summary: ", summary_path)
  print(summary_dt[, .(sample, cluster_id, median_ccf, n_mutations,
                        n_segments_with_muts)])
}

message("\nDone.")
