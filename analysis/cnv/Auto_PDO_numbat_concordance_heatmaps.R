####################
# Auto_PDO_numbat_concordance_heatmaps.R
#
# Compare expression-derived InferCNA subclones against Numbat haplotype-aware
# clone calls. Writes one PDF page per sample with matched-cell CNV heatmaps:
# left = InferCNA-derived profile/cut, right = Numbat profile/cut.
####################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(scales)
})

root_dir <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
out_root <- file.path(root_dir, "PDOs_outs")
setwd(out_root)

out_dir <- "Auto_PDO_numbat"
plot_dir <- file.path(out_dir, "concordance")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(out_dir, "Auto_PDO_numbat_manifest.csv")
infer_outs_path <- "cnv/Auto_PDO_infercna_target_outs_Carroll_2023.rds"
infer_cells_path <- "Auto_PDO_cnv_subclone_mp/Auto_PDO_cnv_subclone_cells.csv"
gene_order_path <- "/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt"
mp_score_path <- "Auto_PDO_mp_adj_noreg.rds"

for (path in c(manifest_path, infer_outs_path, infer_cells_path, gene_order_path)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

max_plot_cells <- 1200L
cna_colour_limit <- 0.15
numbat_colour_limit <- 1
set.seed(123)

####################
# The loop explicitly opens one page per sample with grid.newpage(). Keep
# grid.arrange() on that page rather than letting it open a second blank page.
grid_arrange_original <- gridExtra::grid.arrange
grid.arrange <- function(..., newpage = FALSE) {
  grid_arrange_original(..., newpage = newpage)
}
####################

manifest <- fread(manifest_path)
infer_cna <- readRDS(infer_outs_path)
infer_cells <- fread(infer_cells_path) %>%
  rename(cell_id = .data$cell, infercna_subclone = .data$subclone)

gene_order <- fread(
  gene_order_path,
  header = FALSE,
  col.names = c("gene", "chromosome", "start", "end")
)
####################
# The local hg38 gene-order table can contain occasional extra empty fields.
# Re-read with fill=TRUE and keep the first four columns so chromosome binning
# uses the full genome rather than stopping at the first irregular line.
gene_order <- fread(gene_order_path, header = FALSE, fill = TRUE)
gene_order <- gene_order[, seq_len(4), with = FALSE]
setnames(gene_order, c("gene", "chromosome", "start", "end"))
gene_order <- gene_order[!is.na(gene) & nzchar(gene)]
####################
chrom_levels <- c(paste0("chr", 1:22), "chrX")
gene_order <- gene_order %>%
  filter(.data$chromosome %in% chrom_levels) %>%
  mutate(start = as.numeric(.data$start), end = as.numeric(.data$end)) %>%
  filter(is.finite(.data$start), is.finite(.data$end)) %>%
  mutate(chromosome = factor(.data$chromosome, levels = chrom_levels)) %>%
  arrange(.data$chromosome, .data$start)

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

centromere_pos <- c(
  chr1 = 121700000, chr2 = 91800000, chr3 = 87900000, chr4 = 50600000,
  chr5 = 48400000, chr6 = 61000000, chr7 = 59900000, chr8 = 45600000,
  chr9 = 49000000, chr10 = 40200000, chr11 = 53400000, chr12 = 35500000,
  chr13 = 17700000, chr14 = 17200000, chr15 = 19000000, chr16 = 36800000,
  chr17 = 25100000, chr18 = 18500000, chr19 = 26200000, chr20 = 28100000,
  chr21 = 12000000, chr22 = 15000000, chrX = 61000000
)

make_palette <- function(values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(character(0))
  base <- suppressWarnings(brewer.pal(max(3, min(12, length(values))), palette))
  setNames(colorRampPalette(base)(length(values)), values)
}

complete_palette <- function(cols, values, palette = "Set3") {
  values <- sort(unique(as.character(values)))
  values <- values[!is.na(values) & nzchar(values)]
  cols <- cols[!is.na(names(cols)) & nzchar(names(cols))]
  missing_values <- setdiff(values, names(cols))
  if (length(missing_values) > 0) cols <- c(cols, make_palette(missing_values, palette))
  cols[values]
}

comb2 <- function(x) x * (x - 1) / 2

adjusted_rand <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- as.character(x[ok])
  y <- as.character(y[ok])
  n <- length(x)
  if (n < 2) return(NA_real_)
  tab <- table(x, y)
  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- sum_i * sum_j / total
  max_index <- (sum_i + sum_j) / 2
  denom <- max_index - expected
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (sum_ij - expected) / denom
}

normalised_mi <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- as.character(x[ok])
  y <- as.character(y[ok])
  n <- length(x)
  if (n < 2) return(NA_real_)
  tab <- table(x, y)
  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- pxy > 0
  mi <- sum(pxy[nz] * log(pxy[nz] / outer(px, py)[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx <= 0 || hy <= 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

cluster_purity <- function(reference, query) {
  ok <- !is.na(reference) & !is.na(query)
  tab <- table(as.character(query[ok]), as.character(reference[ok]))
  if (sum(tab) == 0) return(NA_real_)
  sum(apply(tab, 1, max)) / sum(tab)
}

read_numbat_clone_post <- function(row) {
  sample_id <- row$sample
  clone_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_clone_post.csv"))
  map_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_cell_map.csv"))
  if (!file.exists(clone_file) || !file.exists(map_file)) return(NULL)
  clone_post <- fread(clone_file)
  cell_map <- fread(map_file)
  if (!("cell" %in% colnames(clone_post)) || !("clone_opt" %in% colnames(clone_post))) return(NULL)
  clone_post <- clone_post %>%
    left_join(cell_map %>% select(.data$numbat_cell, .data$cell_id, .data$raw_barcode), by = c("cell" = "numbat_cell"))
  missing_cell <- is.na(clone_post$cell_id)
  if (any(missing_cell)) {
    raw_guess <- sub(paste0("^", sample_id, "_"), "", clone_post$cell[missing_cell])
    fallback <- cell_map$cell_id[match(raw_guess, cell_map$raw_barcode)]
    clone_post$cell_id[missing_cell] <- fallback
  }
  clone_post$numbat_p_cnv <- if ("p_cnv" %in% colnames(clone_post)) clone_post$p_cnv else NA_real_
  clone_post$numbat_p_cnv_expr <- if ("p_cnv_x" %in% colnames(clone_post)) clone_post$p_cnv_x else NA_real_
  clone_post$numbat_p_cnv_allele <- if ("p_cnv_y" %in% colnames(clone_post)) clone_post$p_cnv_y else NA_real_
  clone_post %>%
    transmute(
      cell_id = .data$cell_id,
      sample = sample_id,
      numbat_cell = .data$cell,
      numbat_clone = paste0("Numbat clone ", .data$clone_opt),
      numbat_clone_opt = as.character(.data$clone_opt),
      numbat_p_cnv = .data$numbat_p_cnv,
      numbat_p_cnv_expr = .data$numbat_p_cnv_expr,
      numbat_p_cnv_allele = .data$numbat_p_cnv_allele
    ) %>%
    filter(!is.na(.data$cell_id))
}

score_numbat_joint <- function(joint_post) {
  p_cols <- colnames(joint_post)
  if (all(c("p_amp", "p_bamp", "p_del", "p_bdel") %in% p_cols)) {
    return(
      as.numeric(joint_post$p_amp) +
        as.numeric(joint_post$p_bamp) -
        as.numeric(joint_post$p_del) -
        as.numeric(joint_post$p_bdel)
    )
  }
  state <- as.character(joint_post$cnv_state)
  base <- case_when(
    state %in% c("amp", "bamp") ~ 1,
    state %in% c("del", "bdel") ~ -1,
    state %in% c("loh") ~ 0.5,
    TRUE ~ 0
  )
  prob <- if ("p_cnv" %in% p_cols) as.numeric(joint_post$p_cnv) else 1
  base * prob
}

read_numbat_matrix <- function(row, cells_use) {
  sample_id <- row$sample
  joint_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_joint_post.csv.gz"))
  map_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_cell_map.csv"))
  if (!file.exists(joint_file) || !file.exists(map_file)) return(NULL)
  joint_post <- fread(joint_file)
  cell_map <- fread(map_file)
  joint_post <- joint_post %>%
    left_join(cell_map %>% select(.data$numbat_cell, .data$cell_id), by = c("cell" = "numbat_cell")) %>%
    filter(.data$cell_id %in% cells_use)
  if (nrow(joint_post) == 0) return(NULL)
  joint_post$score <- score_numbat_joint(joint_post)
  joint_post$CHROM <- paste0("chr", sub("^chr", "", as.character(joint_post$CHROM)))
  joint_post$CHROM <- factor(joint_post$CHROM, levels = chrom_levels)
  if (!("seg_start" %in% colnames(joint_post))) joint_post$seg_start <- seq_len(nrow(joint_post))
  if (!("seg_end" %in% colnames(joint_post))) joint_post$seg_end <- joint_post$seg_start
  joint_post$row_key <- paste(joint_post$CHROM, joint_post$seg, joint_post$seg_start, joint_post$seg_end, sep = "_")
  joint_post <- joint_post %>%
    filter(!is.na(.data$CHROM)) %>%
    group_by(.data$row_key, .data$CHROM, .data$seg_start, .data$seg_end, .data$cell_id) %>%
    summarise(score = mean(.data$score, na.rm = TRUE), .groups = "drop")
  wide <- dcast(
    as.data.table(joint_post),
    row_key + CHROM + seg_start + seg_end ~ cell_id,
    value.var = "score",
    fill = 0
  )
  wide <- wide[order(CHROM, seg_start, seg_end)]
  present_cells <- intersect(cells_use, colnames(wide))
  if (length(present_cells) == 0) return(NULL)
  mat <- as.matrix(wide[, present_cells, with = FALSE])
  rownames(mat) <- wide$row_key
  list(mat = mat, chr = as.character(wide$CHROM))
}

make_infer_binned <- function(cells_use, bin_size = 100L) {
  cells_use <- intersect(cells_use, colnames(infer_cna))
  common_genes <- intersect(rownames(infer_cna), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE]
  go <- go[order(go$chromosome, go$start), , drop = FALSE]
  mat <- as.matrix(infer_cna[go$gene, cells_use, drop = FALSE])
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]
  go$.row <- seq_len(nrow(go))
  go <- go %>%
    group_by(.data$chromosome) %>%
    mutate(bin = paste0(.data$chromosome, "_", ((row_number() - 1L) %/% bin_size) + 1L)) %>%
    ungroup()
  bins <- split(go$.row, factor(go$bin, levels = unique(go$bin)))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  list(mat = binned, chr = sub("_.*$", "", rownames(binned)))
}

####################
# InferCNA output columns use sample__barcode, while downstream state/subclone
# tables use sample_barcode. Resolve the matrix columns, then restore the
# downstream cell IDs so InferCNA and Numbat panels share a matched-cell axis.
resolve_infer_columns <- function(cells_use) {
  infer_cols <- colnames(infer_cna)
  sample_ids <- manifest$sample
  resolved <- vapply(cells_use, function(cell_id) {
    candidates <- cell_id
    for (sample_id in sample_ids) {
      prefix <- paste0(sample_id, "_")
      if (startsWith(cell_id, prefix)) {
        raw_barcode <- sub(paste0("^", prefix), "", cell_id)
        candidates <- c(
          paste(sample_id, raw_barcode, sep = "__"),
          cell_id,
          raw_barcode
        )
        break
      }
    }
    hit <- candidates[candidates %in% infer_cols]
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }, character(1))
  resolved
}

make_infer_binned <- function(cells_use, bin_size = 100L) {
  infer_lookup <- resolve_infer_columns(cells_use)
  keep_cells <- !is.na(infer_lookup)
  cells_use <- cells_use[keep_cells]
  infer_lookup <- infer_lookup[keep_cells]
  if (length(cells_use) == 0) stop("No InferCNA matrix columns matched requested cells.")

  common_genes <- intersect(rownames(infer_cna), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE]
  go <- go[order(go$chromosome, go$start), , drop = FALSE]
  mat <- as.matrix(infer_cna[go$gene, infer_lookup, drop = FALSE])
  colnames(mat) <- cells_use
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]
  go$.row <- seq_len(nrow(go))
  go <- go %>%
    group_by(.data$chromosome) %>%
    mutate(bin = paste0(.data$chromosome, "_", ((row_number() - 1L) %/% bin_size) + 1L)) %>%
    ungroup()
  bins <- split(go$.row, factor(go$bin, levels = unique(go$bin)))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  colnames(binned) <- cells_use
  list(mat = binned, chr = sub("_.*$", "", rownames(binned)))
}
####################

####################
# Use a shared genome-wide gene-bin scaffold for both InferCNA and Numbat. Numbat
# only emits posterior rows for inferred CNV segments; bins with no emitted
# segment are shown as neutral zero so all chromosomes remain visible.
make_genome_bins <- function(bin_size = 100L) {
  go <- gene_order %>%
    mutate(.row = row_number()) %>%
    group_by(.data$chromosome) %>%
    mutate(bin_index = ((row_number() - 1L) %/% bin_size) + 1L,
           bin = paste0(.data$chromosome, "_", .data$bin_index)) %>%
    ungroup()
  bins <- go %>%
    group_by(.data$bin, .data$chromosome, .data$bin_index) %>%
    summarise(start = min(.data$start, na.rm = TRUE),
              end = max(.data$end, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(.data$chromosome, .data$bin_index)
  bins
}

genome_bins <- make_genome_bins()
####################
# Chromosome-arm labels are used for the right-side scatter plot. Use the same
# fixed centromere positions as the InferCNA subclone workflow so arm means are
# computed consistently across the two CNA sources.
genome_bins <- genome_bins %>%
  mutate(
    midpoint = (.data$start + .data$end) / 2,
    arm = ifelse(.data$midpoint <= centromere_pos[as.character(.data$chromosome)], "p", "q"),
    arm_label = paste0(as.character(.data$chromosome), .data$arm)
  )
####################

make_infer_binned <- function(cells_use) {
  infer_lookup <- resolve_infer_columns(cells_use)
  keep_cells <- !is.na(infer_lookup)
  cells_use <- cells_use[keep_cells]
  infer_lookup <- infer_lookup[keep_cells]
  if (length(cells_use) == 0) stop("No InferCNA matrix columns matched requested cells.")

  common_genes <- intersect(rownames(infer_cna), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE] %>%
    mutate(.row = row_number()) %>%
    group_by(.data$chromosome) %>%
    mutate(bin_index = ((row_number() - 1L) %/% 100L) + 1L,
           bin = paste0(.data$chromosome, "_", .data$bin_index)) %>%
    ungroup()
  mat <- as.matrix(infer_cna[go$gene, infer_lookup, drop = FALSE])
  colnames(mat) <- cells_use
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]

  binned <- matrix(0, nrow = nrow(genome_bins), ncol = length(cells_use),
                   dimnames = list(genome_bins$bin, cells_use))
  bins <- split(seq_len(nrow(go)), factor(go$bin, levels = genome_bins$bin))
  bins <- bins[lengths(bins) > 0]
  for (bin_name in names(bins)) {
    binned[bin_name, ] <- colMeans(mat[bins[[bin_name]], , drop = FALSE], na.rm = TRUE)
  }
  list(mat = binned, chr = as.character(genome_bins$chromosome))
}

####################
# Exact InferCNA preparation from Auto_PDO_cnv_subclone_mp_heatmap.R: keep the
# finite expression-CNA matrix, filter out the lowest third of genes by CNA
# signal, then bin the filtered profile for plotting.
prepare_infer_cna_matrix <- function(cells_use) {
  infer_lookup <- resolve_infer_columns(cells_use)
  keep_cells <- !is.na(infer_lookup)
  cells_use <- cells_use[keep_cells]
  infer_lookup <- infer_lookup[keep_cells]
  if (length(cells_use) == 0) stop("No InferCNA matrix columns matched requested cells.")
  common_genes <- intersect(rownames(infer_cna), gene_order$gene)
  go <- gene_order[match(common_genes, gene_order$gene), , drop = FALSE]
  go <- go[order(go$chromosome, go$start), , drop = FALSE]
  mat <- as.matrix(infer_cna[go$gene, infer_lookup, drop = FALSE])
  colnames(mat) <- cells_use
  keep <- rowSums(is.finite(mat)) == ncol(mat)
  mat <- mat[keep, , drop = FALSE]
  go <- go[keep, , drop = FALSE]
  signal <- rowMeans(abs(mat), na.rm = TRUE)
  keep_signal <- signal >= as.numeric(quantile(signal, probs = 1 / 3, na.rm = TRUE))
  list(mat = mat[keep_signal, , drop = FALSE], gene_order = go[keep_signal, , drop = FALSE])
}

make_binned_cna <- function(cna_mat, go, bin_size = 100L) {
  go2 <- go %>%
    mutate(.row = seq_len(n())) %>%
    group_by(.data$chromosome) %>%
    mutate(bin_index = ((row_number() - 1L) %/% bin_size) + 1L,
           bin = paste0(.data$chromosome, "_", .data$bin_index)) %>%
    ungroup()
  bin_levels <- unique(go2$bin)
  bins <- split(go2$.row, factor(go2$bin, levels = bin_levels))
  binned <- do.call(rbind, lapply(bins, function(ix) colMeans(cna_mat[ix, , drop = FALSE], na.rm = TRUE)))
  rownames(binned) <- names(bins)
  bins_df <- go2 %>%
    group_by(.data$bin, .data$chromosome, .data$bin_index) %>%
    summarise(start = min(.data$start, na.rm = TRUE),
              end = max(.data$end, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(
      midpoint = (.data$start + .data$end) / 2,
      arm = ifelse(.data$midpoint <= centromere_pos[as.character(.data$chromosome)], "p", "q"),
      arm_label = paste0(as.character(.data$chromosome), .data$arm)
    )
  bins_df <- bins_df[match(rownames(binned), bins_df$bin), , drop = FALSE]
  list(mat = binned, chr = sub("_.*$", "", rownames(binned)), bins = bins_df)
}
####################

read_numbat_matrix <- function(row, cells_use, bins_df = genome_bins) {
  sample_id <- row$sample
  joint_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_joint_post.csv.gz"))
  map_file <- file.path(row$numbat_dir, paste0("Auto_", sample_id, "_numbat_cell_map.csv"))
  if (!file.exists(joint_file) || !file.exists(map_file)) return(NULL)
  joint_post <- fread(joint_file)
  cell_map <- fread(map_file)

  available_cells <- intersect(cells_use, cell_map$cell_id)
  if (length(available_cells) == 0) return(NULL)
  mat <- matrix(0, nrow = nrow(bins_df), ncol = length(available_cells),
                dimnames = list(bins_df$bin, available_cells))

  joint_post <- joint_post %>%
    left_join(cell_map %>% select(.data$numbat_cell, .data$cell_id), by = c("cell" = "numbat_cell")) %>%
    filter(.data$cell_id %in% available_cells)
  if (nrow(joint_post) == 0) {
    return(list(mat = mat, chr = as.character(bins_df$chromosome), bins = bins_df))
  }

  joint_post$score <- score_numbat_joint(joint_post)
  joint_post$CHROM <- paste0("chr", sub("^chr", "", as.character(joint_post$CHROM)))
  joint_post$seg_start <- as.numeric(joint_post$seg_start)
  joint_post$seg_end <- as.numeric(joint_post$seg_end)
  joint_post <- joint_post %>%
    filter(.data$CHROM %in% chrom_levels,
           is.finite(.data$seg_start),
           is.finite(.data$seg_end),
           is.finite(.data$score)) %>%
    group_by(.data$CHROM, .data$seg, .data$seg_start, .data$seg_end, .data$cell_id) %>%
    summarise(score = mean(.data$score, na.rm = TRUE), .groups = "drop")
  if (nrow(joint_post) == 0) {
    return(list(mat = mat, chr = as.character(bins_df$chromosome), bins = bins_df))
  }

  segment_keys <- joint_post %>%
    distinct(.data$CHROM, .data$seg, .data$seg_start, .data$seg_end)
  for (i in seq_len(nrow(segment_keys))) {
    seg_row <- segment_keys[i, ]
    bin_ix <- which(
      as.character(bins_df$chromosome) == seg_row$CHROM &
        bins_df$end >= seg_row$seg_start &
        bins_df$start <= seg_row$seg_end
    )
    if (length(bin_ix) == 0) next
    seg_scores <- joint_post %>%
      filter(.data$CHROM == seg_row$CHROM,
             .data$seg == seg_row$seg,
             .data$seg_start == seg_row$seg_start,
             .data$seg_end == seg_row$seg_end)
    cell_ix <- match(seg_scores$cell_id, colnames(mat))
    cell_keep <- is.finite(cell_ix)
    if (!any(cell_keep)) next
    mat[bin_ix, cell_ix[cell_keep]] <- matrix(
      rep(seg_scores$score[cell_keep], each = length(bin_ix)),
      nrow = length(bin_ix)
    )
  }

  list(mat = mat, chr = as.character(bins_df$chromosome), bins = bins_df)
}
####################

order_cells_by_group <- function(mat, group) {
  group <- as.character(group)
  names(group) <- colnames(mat)
  split_cells <- split(colnames(mat), factor(group, levels = unique(group)))
  unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= 2) return(cells)
    d <- dist(t(mat[, cells, drop = FALSE]))
    hc <- hclust(d, method = "ward.D2")
    cells[hc$order]
  }), use.names = FALSE)
}

####################
# Some matrix construction routes can drop column names when only one row/column
# survives a subset. Restore/validate names before within-clone clustering.
order_cells_by_group <- function(mat, group) {
  if (is.null(dim(mat))) mat <- matrix(mat, nrow = 1)
  group <- as.character(group)
  if (is.null(names(group)) || any(is.na(names(group))) || any(!nzchar(names(group)))) {
    names(group) <- colnames(mat)
  }
  if (is.null(colnames(mat)) || any(is.na(colnames(mat))) || any(!nzchar(colnames(mat)))) {
    if (length(group) == ncol(mat) && !is.null(names(group))) {
      colnames(mat) <- names(group)
    } else {
      stop("Heatmap matrix is missing column names and cannot be matched to metadata.")
    }
  }
  group <- group[colnames(mat)]
  group[is.na(group) | !nzchar(group)] <- "NA"
  split_cells <- split(colnames(mat), factor(group, levels = unique(group)))
  unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= 2) return(cells)
    mat_sub <- mat[, cells, drop = FALSE]
    if (any(!is.finite(mat_sub))) mat_sub[!is.finite(mat_sub)] <- 0
    if (nrow(mat_sub) < 2 || ncol(mat_sub) < 3) return(cells)
    d <- dist(t(mat_sub))
    if (length(d) == 0 || any(!is.finite(d))) return(cells)
    hc <- hclust(d, method = "ward.D2")
    cells[hc$order]
  }), use.names = FALSE)
}
####################

sample_cells_for_plot <- function(cells, infer_group, numbat_group) {
  if (length(cells) <= max_plot_cells) return(cells)
  split_cells <- split(cells, factor(infer_group[cells], levels = unique(infer_group[cells])))
  target <- pmax(20L, floor(max_plot_cells * lengths(split_cells) / length(cells)))
  target <- pmin(target, lengths(split_cells))
  out <- unlist(mapply(function(x, n) sample(x, n), split_cells, target, SIMPLIFY = FALSE), use.names = FALSE)
  if (length(out) > max_plot_cells) out <- sample(out, max_plot_cells)
  out
}

make_heatmap_grob <- function(mat, row_chr, meta_plot, group_col, heatmap_name, title, colour_limit,
                              palette = "Set2", show_annotation_legend = TRUE, cell_order = NULL) {
  cells <- rownames(meta_plot)
  mat <- mat[, cells, drop = FALSE]
  group <- meta_plot[[group_col]]
  order <- if (!is.null(cell_order)) {
    intersect(cell_order, colnames(mat))
  } else {
    order_cells_by_group(mat, group)
  }
  mat <- mat[, order, drop = FALSE]
  meta_plot <- meta_plot[order, , drop = FALSE]
  group <- meta_plot[[group_col]]

  lim <- colour_limit
  group_cols <- complete_palette(make_palette(group, palette), group, palette)
  infer_cols <- complete_palette(make_palette(meta_plot$infercna_subclone, "Set2"), meta_plot$infercna_subclone, "Set2")
  numbat_cols <- complete_palette(make_palette(meta_plot$numbat_clone, "Dark2"), meta_plot$numbat_clone, "Dark2")
  local_state_cols <- complete_palette(state_cols, meta_plot$state_label, "Set3")
  local_mp_cols <- complete_palette(mp_cols, meta_plot$top_mp_label, "Paired")

  top_ha <- HeatmapAnnotation(
    InferCNA = meta_plot$infercna_subclone,
    Numbat = meta_plot$numbat_clone,
    State = meta_plot$state_label,
    TopMP = meta_plot$top_mp_label,
    col = list(
      InferCNA = infer_cols,
      Numbat = numbat_cols,
      State = local_state_cols,
      TopMP = local_mp_cols
    ),
    annotation_name_side = "left",
    show_annotation_name = TRUE,
    simple_anno_size = unit(4, "mm"),
    show_legend = TRUE,
    na_col = "grey90"
  )

  row_chr_factor <- factor(row_chr, levels = unique(row_chr))
  chr_cols <- setNames(rep(c("#E6E6E6", "#BDBDBD"), length.out = length(levels(row_chr_factor))), levels(row_chr_factor))
  left_ha <- rowAnnotation(
    Chr = row_chr_factor,
    col = list(Chr = chr_cols),
    show_annotation_name = FALSE,
    show_legend = FALSE,
    width = unit(3, "mm")
  )

  ht <- Heatmap(
    mat,
    name = heatmap_name,
    col = colorRamp2(c(-lim, 0, lim), c("#2166AC", "white", "#B2182B")),
    left_annotation = left_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = factor(group, levels = unique(group)),
    row_split = row_chr_factor,
    row_gap = unit(0, "mm"),
    column_gap = unit(1.2, "mm"),
    show_column_names = FALSE,
    show_row_names = FALSE,
    show_column_dend = FALSE,
    top_annotation = top_ha,
    column_title = title,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_gp = gpar(fontsize = 7),
    use_raster = TRUE,
    raster_quality = 4,
    border = FALSE,
    rect_gp = gpar(col = NA),
    heatmap_legend_param = list(title = heatmap_name, legend_height = unit(26, "mm"))
  )

  grid.grabExpr(draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    show_annotation_legend = show_annotation_legend
  ))
}

make_matching_bar_plot <- function(meta_all) {
  plot_df <- meta_all %>%
    count(.data$infercna_subclone, .data$numbat_clone, name = "n") %>%
    group_by(.data$infercna_subclone) %>%
    mutate(frac = .data$n / sum(.data$n)) %>%
    ungroup()
  numbat_cols <- complete_palette(make_palette(plot_df$numbat_clone, "Dark2"), plot_df$numbat_clone, "Dark2")
  ggplot(plot_df, aes(x = .data$infercna_subclone, y = .data$frac, fill = .data$numbat_clone)) +
    geom_col(color = "white", linewidth = 0.2, width = 0.78) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = numbat_cols, drop = FALSE) +
    labs(title = "Subclone overlap", x = "InferCNA subclone", y = "Numbat fraction", fill = "Numbat") +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
      legend.position = "bottom",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 5),
      legend.key.size = unit(2.5, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

arm_mean_matrix <- function(mat, bins_df) {
  stopifnot(nrow(mat) == nrow(bins_df))
  arm_levels <- unique(bins_df$arm_label)
  out <- do.call(rbind, lapply(arm_levels, function(arm) {
    rows <- which(bins_df$arm_label == arm)
    colMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(out) <- arm_levels
  out[!is.finite(out)] <- 0
  out
}

make_arm_scatter_plot <- function(infer_mat, numbat_mat, bins_df) {
  common_cells <- intersect(colnames(infer_mat), colnames(numbat_mat))
  infer_arm <- arm_mean_matrix(infer_mat[, common_cells, drop = FALSE], bins_df)
  numbat_arm <- arm_mean_matrix(numbat_mat[, common_cells, drop = FALSE], bins_df)
  common_arms <- intersect(rownames(infer_arm), rownames(numbat_arm))
  scatter_df <- data.frame(
    cell_id = rep(common_cells, each = length(common_arms)),
    arm = rep(common_arms, times = length(common_cells)),
    infercna = as.vector(infer_arm[common_arms, common_cells, drop = FALSE]),
    numbat = as.vector(numbat_arm[common_arms, common_cells, drop = FALSE]),
    stringsAsFactors = FALSE
  )
  scatter_df <- scatter_df[is.finite(scatter_df$infercna) & is.finite(scatter_df$numbat), , drop = FALSE]
  rho <- suppressWarnings(cor(scatter_df$infercna, scatter_df$numbat, method = "spearman", use = "complete.obs"))
  subtitle <- paste0("chromosome-arm means, cells x arms=", nrow(scatter_df), ", Spearman rho=", signif(rho, 3))
  x_lim <- as.numeric(quantile(scatter_df$infercna, c(0.01, 0.99), na.rm = TRUE))
  y_lim <- as.numeric(quantile(scatter_df$numbat, c(0.01, 0.99), na.rm = TRUE))
  if (!all(is.finite(x_lim)) || diff(x_lim) <= 0) x_lim <- range(scatter_df$infercna, finite = TRUE)
  if (!all(is.finite(y_lim)) || diff(y_lim) <= 0) y_lim <- range(scatter_df$numbat, finite = TRUE)
  ggplot(scatter_df, aes(.data$infercna, .data$numbat)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.25) +
    geom_point(color = "#2B8CBE", alpha = 0.12, size = 0.18) +
    geom_smooth(method = "lm", se = FALSE, color = "#B2182B", linewidth = 0.35) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(title = "Arm CNA agreement", subtitle = subtitle, x = "InferCNA arm mean", y = "Numbat arm mean") +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 6),
      plot.margin = margin(2, 2, 2, 2)
    )
}

clone_posts <- lapply(seq_len(nrow(manifest)), function(i) read_numbat_clone_post(manifest[i]))
clone_posts <- bind_rows(clone_posts)
if (nrow(clone_posts) == 0) stop("No Numbat clone summaries were found.")

cell_df <- infer_cells %>%
  inner_join(clone_posts, by = c("cell_id", "sample")) %>%
  mutate(
    infercna_subclone = as.character(.data$infercna_subclone),
    numbat_clone = as.character(.data$numbat_clone),
    top_mp_label = ifelse(is.na(.data$top_mp_label), "NA", as.character(.data$top_mp_label)),
    state_label = ifelse(is.na(.data$state_label), "NA", as.character(.data$state_label))
  )

fwrite(cell_df, file.path(plot_dir, "Auto_PDO_numbat_cell_concordance.csv"))

summary_rows <- cell_df %>%
  group_by(.data$sample) %>%
  summarise(
    n_common_cells = n(),
    n_infercna_subclones = n_distinct(.data$infercna_subclone),
    n_numbat_clones = n_distinct(.data$numbat_clone),
    adjusted_rand = adjusted_rand(.data$infercna_subclone, .data$numbat_clone),
    normalised_mi = normalised_mi(.data$infercna_subclone, .data$numbat_clone),
    purity_numbat_given_infercna = cluster_purity(.data$infercna_subclone, .data$numbat_clone),
    purity_infercna_given_numbat = cluster_purity(.data$numbat_clone, .data$infercna_subclone),
    median_numbat_p_cnv = median(.data$numbat_p_cnv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(.data$sample)
fwrite(summary_rows, file.path(plot_dir, "Auto_PDO_numbat_infercna_concordance_summary.csv"))

contingency <- cell_df %>%
  count(.data$sample, .data$infercna_subclone, .data$numbat_clone, name = "n_cells") %>%
  group_by(.data$sample, .data$infercna_subclone) %>%
  mutate(frac_of_infercna = .data$n_cells / sum(.data$n_cells)) %>%
  ungroup() %>%
  group_by(.data$sample, .data$numbat_clone) %>%
  mutate(frac_of_numbat = .data$n_cells / sum(.data$n_cells)) %>%
  ungroup()
fwrite(contingency, file.path(plot_dir, "Auto_PDO_numbat_infercna_contingency.csv"))

state_summary <- cell_df %>%
  count(.data$sample, .data$numbat_clone, .data$state_label, name = "n_cells") %>%
  group_by(.data$sample, .data$numbat_clone) %>%
  mutate(frac = .data$n_cells / sum(.data$n_cells)) %>%
  ungroup()
fwrite(state_summary, file.path(plot_dir, "Auto_PDO_numbat_clone_state_summary.csv"))

topmp_summary <- cell_df %>%
  count(.data$sample, .data$numbat_clone, .data$top_mp, .data$top_mp_label, name = "n_cells") %>%
  group_by(.data$sample, .data$numbat_clone) %>%
  mutate(frac = .data$n_cells / sum(.data$n_cells)) %>%
  ungroup()
fwrite(topmp_summary, file.path(plot_dir, "Auto_PDO_numbat_clone_topmp_summary.csv"))

if (file.exists(mp_score_path)) {
  mp_scores <- as.matrix(readRDS(mp_score_path))
  score_cells <- intersect(cell_df$cell_id, rownames(mp_scores))
  if (length(score_cells) > 0) {
    mp_base <- cell_df %>%
      filter(.data$cell_id %in% score_cells) %>%
      select(.data$cell_id, .data$sample, .data$numbat_clone)
    mp_summary_input <- bind_cols(
      mp_base,
      as.data.frame(mp_scores[mp_base$cell_id, , drop = FALSE])
    )
    mp_summary <- mp_summary_input %>%
      group_by(.data$sample, .data$numbat_clone) %>%
      summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
    fwrite(mp_summary, file.path(plot_dir, "Auto_PDO_numbat_clone_mp_score_summary.csv"))
  }
}

output_pdf <- file.path(plot_dir, "Auto_PDO_numbat_infercna_matched_heatmaps.pdf")
pdf(output_pdf, width = 21, height = 11.5, useDingbats = FALSE)
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  sample_id <- row$sample
  sample_cells <- cell_df %>% filter(.data$sample == sample_id)
  if (nrow(sample_cells) < 20) {
    grid.newpage()
    grid.text(paste0(sample_id, "\nToo few matched InferCNA/Numbat cells for heatmap."), gp = gpar(fontsize = 16))
    next
  }

  infer_group <- setNames(sample_cells$infercna_subclone, sample_cells$cell_id)
  numbat_group <- setNames(sample_cells$numbat_clone, sample_cells$cell_id)
  infer_prepped_all <- prepare_infer_cna_matrix(sample_cells$cell_id)
  infer_binned_all <- make_binned_cna(infer_prepped_all$mat, infer_prepped_all$gene_order)
  numbat_binned_all <- read_numbat_matrix(row, colnames(infer_binned_all$mat), infer_binned_all$bins)
  if (is.null(numbat_binned_all) || ncol(numbat_binned_all$mat) < 20 || nrow(numbat_binned_all$mat) < 2) {
    grid.newpage()
    grid.text(paste0(sample_id, "\nMissing usable Numbat joint posterior matrix."), gp = gpar(fontsize = 16))
    next
  }

  common_all_cells <- Reduce(intersect, list(
    sample_cells$cell_id,
    colnames(infer_binned_all$mat),
    colnames(numbat_binned_all$mat)
  ))
  sample_cells <- sample_cells[match(common_all_cells, sample_cells$cell_id), ]
  infer_group <- setNames(sample_cells$infercna_subclone, sample_cells$cell_id)
  numbat_group <- setNames(sample_cells$numbat_clone, sample_cells$cell_id)

  set.seed(42)
  plot_cells <- sample_cells_for_plot(sample_cells$cell_id, infer_group, numbat_group)
  infer_order <- order_cells_by_group(
    infer_prepped_all$mat[, plot_cells, drop = FALSE],
    factor(infer_group[plot_cells], levels = unique(infer_group[plot_cells]))
  )
  sample_cells <- sample_cells[match(plot_cells, sample_cells$cell_id), ]
  meta_plot <- as.data.frame(sample_cells)
  rownames(meta_plot) <- meta_plot$cell_id

  infer_binned <- list(
    mat = infer_binned_all$mat[, plot_cells, drop = FALSE],
    chr = infer_binned_all$chr
  )
  numbat_binned <- list(
    mat = numbat_binned_all$mat[, plot_cells, drop = FALSE],
    chr = numbat_binned_all$chr
  )

  ####################
  # Preserve matched-cell names before intersecting panels; a zero-name matrix
  # would otherwise make the plotting step fail without identifying the sample.
  if (is.null(colnames(infer_binned$mat)) && ncol(infer_binned$mat) == length(plot_cells)) {
    colnames(infer_binned$mat) <- plot_cells
  }
  if (is.null(colnames(numbat_binned$mat)) && ncol(numbat_binned$mat) == length(plot_cells)) {
    colnames(numbat_binned$mat) <- plot_cells
  }
  if (is.null(colnames(infer_binned$mat)) || is.null(colnames(numbat_binned$mat))) {
    stop("Missing heatmap column names for sample: ", sample_id)
  }
  message(
    "Rendering ", sample_id,
    ": infer rows=", nrow(infer_binned$mat),
    ", numbat rows=", nrow(numbat_binned$mat),
    ", heatmap cells=", length(plot_cells),
    ", scatter cells=", length(common_all_cells)
  )
  ####################

  plot_cells <- Reduce(intersect, list(colnames(infer_binned$mat), colnames(numbat_binned$mat), rownames(meta_plot)))
  meta_plot <- meta_plot[plot_cells, , drop = FALSE]
  infer_binned$mat <- infer_binned$mat[, plot_cells, drop = FALSE]
  numbat_binned$mat <- numbat_binned$mat[, plot_cells, drop = FALSE]

  stats_row <- summary_rows %>% filter(.data$sample == sample_id)
  page_title <- sprintf(
    "%s | cells=%s | ARI=%.3f | NMI=%.3f",
    sample_id,
    length(common_all_cells),
    stats_row$adjusted_rand[1],
    stats_row$normalised_mi[1]
  )

  left <- make_heatmap_grob(
    infer_binned$mat,
    infer_binned$chr,
    meta_plot,
    "infercna_subclone",
    "CNA",
    "InferCNA subclones",
    cna_colour_limit,
    "Set2"
    ,
    TRUE,
    infer_order
  )
  right <- make_heatmap_grob(
    numbat_binned$mat,
    numbat_binned$chr,
    meta_plot,
    "numbat_clone",
    "CNV posterior",
    "Numbat subclones",
    numbat_colour_limit,
    "Dark2",
    FALSE
  )
  overlap_plot <- make_matching_bar_plot(cell_df %>% filter(.data$sample == sample_id, .data$cell_id %in% common_all_cells))
  scatter_plot <- make_arm_scatter_plot(
    infer_binned_all$mat[, common_all_cells, drop = FALSE],
    numbat_binned_all$mat[, common_all_cells, drop = FALSE],
    infer_binned_all$bins
  )

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    nrow = 3,
    ncol = 3,
    heights = unit.c(unit(0.35, "in"), unit(1, "null"), unit(1, "null")),
    widths = unit(c(4.15, 4.15, 2.25), "null")
  )))
  grid.text(page_title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:3),
            gp = gpar(fontsize = 13, fontface = "bold"))
  pushViewport(viewport(layout.pos.row = 2:3, layout.pos.col = 1))
  grid.draw(left)
  popViewport()
  pushViewport(viewport(layout.pos.row = 2:3, layout.pos.col = 2))
  grid.draw(right)
  popViewport()
  print(overlap_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
  print(scatter_plot, vp = viewport(layout.pos.row = 3, layout.pos.col = 3))
  popViewport()
}
dev.off()

message("Wrote: ", file.path(out_root, output_pdf))
