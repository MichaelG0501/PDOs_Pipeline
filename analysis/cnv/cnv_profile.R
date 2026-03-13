library(infercna)
library(dplyr)
library(ggplot2)
library(Seurat)

data1 <- readRDS("../matched/patient_H.rds")
data1@meta.data$ident <- rep("patient_H", ncol(data1))
data1_sub <- subset(data1, celltypist == "Ductal")
meta1 <- data1@meta.data[, c("orig.ident", "celltypist", "ident")]
outs1 <- readRDS("../matched/patient_H_outs.rds")
outs1 <- outs1[, rownames(meta1)]

data2 <- tmdata_annotated$SUR680T3_PDO
data2@meta.data$ident <- rep("SUR680", ncol(data2))
data2_sub <- subset(data2, celltypist == "Ductal")
meta2 <- data2@meta.data[, c("orig.ident", "celltypist", "ident")]
outs2 <- readRDS("PDO680_outs.rds")
outs2 <- outs2[, rownames(meta2)]

colnames(outs1) <- paste0(meta1$ident, "_", colnames(outs1))
rownames(meta1) <- colnames(outs1)
colnames(outs2) <- paste0(meta2$ident, "_", colnames(outs2))
rownames(meta2) <- colnames(outs2)
common_rows <- intersect(rownames(outs1), rownames(outs2))
outs1_common <- outs1[common_rows, ]
outs2_common <- outs2[common_rows, ]
outs <- cbind(outs1_common, outs2_common)
meta <- rbind(meta1, meta2)

library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

gene_order <- read.table("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/infercna/hg38_gencode_v27.txt", header=FALSE)
colnames(gene_order) <- c("gene_id", "chromosome", "start", "end")

# Keep only genes in both datasets
common_genes <- intersect(rownames(outs), gene_order$gene_id)
outs <- outs[common_genes, ]
coord <- cnaScatterPlot(outs)
meta$cna_signal <- coord$cna.signal
meta$cna_cor <- coord$cna.cor


outs <- outs[, rownames(meta)[meta$celltypist == "Ductal"]]

#outs <- outs[, cell_summary$cell_id[cell_summary$classification == "malignant"]]
gene_order_filtered <- gene_order %>% filter(gene_id %in% common_genes)

chrom_levels <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
gene_order_filtered$chromosome <- factor(gene_order_filtered$chromosome, levels = chrom_levels)

# Order genes by chromosome and start position
gene_order_sorted <- gene_order_filtered %>%
  arrange(chromosome, start)

outs <- outs[gene_order_sorted$gene_id, ]

bin_size <- 100
n_bins <- ceiling(nrow(outs) / bin_size)

# Create gene bins
binned_mat <- do.call(rbind, lapply(split(1:nrow(outs), ceiling(seq_along(1:nrow(outs))/bin_size)), function(idx) {
  colMeans(outs[idx, , drop=FALSE])
}))

# Assign rownames for bins
rownames(binned_mat) <- paste0("Bin_", seq_len(n_bins))


meta_filtered <- meta[colnames(binned_mat), , drop=FALSE]

library(RColorBrewer)

idents <- sort(unique(na.omit(meta_filtered$ident)))
n_colors <- length(idents)

# Use a high-contrast palette optimized for small numbers
palette_name <- "Set1"  # Good for ≤8; has vivid tones

ident_colors <- setNames(
  brewer.pal(n_colors, palette_name)[seq_along(idents)],
  idents
)

ha <- HeatmapAnnotation(
  orig.ident = meta_filtered$ident,
  col = list(ident = ident_colors),
  show_annotation_name = TRUE
)

# Step 1: Create bin annotation from sorted gene order
gene_order_sorted$bin <- ceiling(seq_len(nrow(gene_order_sorted)) / bin_size)
bin_annotation <- gene_order_sorted %>%
  group_by(bin) %>%
  summarise(
    chr = first(chromosome),
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

# Step 2: Assign bin labels to heatmap
rownames(binned_mat) <- paste0(
  bin_annotation$chr, ":", bin_annotation$start, "-", bin_annotation$end
)

chr_strip <- bin_annotation$chr
chrom_levels_used <- chrom_levels[chrom_levels %in% unique(chr_strip)]
row_split <- factor(chr_strip, levels = chrom_levels_used)

# Step 4: Custom coloring—chr2 and chr11 in red, others in black
custom_chr_colors <- setNames(
  rep("black", length(chrom_levels_used)),
  chrom_levels_used
)
custom_chr_colors[c("chr2", "chr11")] <- "red"

row_anno <- rowAnnotation(
  Chromosome = row_split,
  col = list(Chromosome = custom_chr_colors),
  show_annotation_name = FALSE,
  annotation_name_side = "top",
  annotation_legend_param = list(
    Chromosome = list(
      title = NULL,
      at = NULL,
      labels = NULL
    )
  ),
  gp = gpar(col = NA),
  width = unit(5, "mm")
)


# Distance matrix and hierarchical clustering on transposed matrix (cells)
cell_dist <- dist(t(binned_mat))  # Cells are columns
cell_clust <- hclust(cell_dist, method = "ward.D2")

# Cut tree into desired number of clusters (adjust k as needed)
k <- 10
cell_clusters <- cutree(cell_clust, k = k)

# Add cluster annotation to metadata
meta_filtered$cluster <- factor(cell_clusters[colnames(binned_mat)])

# Assign colors to clusters
cluster_ids <- sort(unique(meta_filtered$cluster))

coord <- cnaScatterPlot(outs)
meta_filtered$cna_signal <- coord$cna.signal
meta_filtered$cna_cor <- coord$cna.cor

valid_cluster_ids <- sort(unique(na.omit(meta_filtered$cluster)))

n_colors <- max(3, min(length(valid_cluster_ids), 8))
cluster_colors <- setNames(
  RColorBrewer::brewer.pal(n_colors, "Set2")[seq_along(valid_cluster_ids)],
  valid_cluster_ids
)
valid_idents <- sort(unique(na.omit(meta_filtered$ident)))

ident_colors <- ident_colors[valid_idents]
ident_colors <- setNames(ident_colors, valid_idents)

cluster_mean_signal <- aggregate(cna_signal ~ cluster, data = meta_filtered, FUN = mean)
cluster_mean_cor <- aggregate(cna_cor ~ cluster, data = meta_filtered, FUN = mean)

# Create named vectors for easy mapping
mean_signal_map <- setNames(cluster_mean_signal$cna_signal, cluster_mean_signal$cluster)
mean_cor_map <- setNames(cluster_mean_cor$cna_cor, cluster_mean_cor$cluster)

# Map values to all cells in meta_filtered
meta_filtered$mean_signal_per_cluster <- mean_signal_map[as.character(meta_filtered$cluster)]
meta_filtered$mean_cor_per_cluster <- mean_cor_map[as.character(meta_filtered$cluster)]

library(circlize)

max_signal <- max(meta_filtered$mean_signal_per_cluster, na.rm = TRUE)
max_cor <- max(meta_filtered$mean_cor_per_cluster, na.rm = TRUE)

col_list <- list(
  orig.ident = ident_colors,
  cluster = cluster_colors,
  mean_cna_signal = colorRamp2(c(0, max_signal), c("white", "red")),
  mean_cna_cor = colorRamp2(c(0, max_cor), c("white", "blue"))
)

meta_filtered$cna_cor_cat <- cut(
  meta_filtered$cna_cor,
  breaks = c(-Inf, 0.15, 0.25, Inf),
  labels = c("<0.15", "0.15–0.25", ">0.25"),
  right = FALSE
)

# Define colors for cna_cor_cat
cna_cor_cat_colors <- setNames(
  c("#D3D3D3", "#87CEFA", "#4169E1"),
  c("<0.15", "0.15–0.25", ">0.25")
)

# Add this to the annotation color list
col_list$cna_cor_cat <- cna_cor_cat_colors

# Update annotation block
ha <- HeatmapAnnotation(
  orig.ident = meta_filtered$ident,
  cluster = meta_filtered$cluster,
  mean_cna_signal = meta_filtered$mean_signal_per_cluster,
  mean_cna_cor = meta_filtered$mean_cor_per_cluster,
  cna_cor_cat = meta_filtered$cna_cor_cat,
  col = col_list,
  annotation_name_side = "left",
  annotation_legend_param = list(
    cluster = list(title = "Cluster"),
    orig.ident = list(title = "Sample"),
    mean_cna_signal = list(title = "Mean CNA Signal"),
    mean_cna_cor = list(title = "Mean CNA Correlation"),
    cna_cor_cat = list(title = "CNA Correlation Category")
  ),
  annotation_height = unit(c(5, 5, 4, 4, 4), "mm"),
  simple_anno_size = unit(2, "mm")
)

# pdf(paste0(unique(meta$ident), "_CNV_profile.pdf"), width = 15, height = 8)
# draw(
#   Heatmap(
#     binned_mat,
#     name = "CNV",
#     cluster_rows = FALSE,
#     cluster_columns = FALSE,
#     column_order = cell_clust$order,
#     show_row_names = FALSE,
#     show_column_names = FALSE,
#     top_annotation = ha,
#     column_title = "CNV Heatmap (Binned by Genomic Position)",
#     row_title_rot = 0,
#     use_raster = FALSE
#   )
# )
# dev.off()


library(zoo)

rolling_window <- 50  # full window
half_window <- floor(rolling_window / 2)

smoothed_mat <- matrix(NA, nrow = nrow(binned_mat), ncol = ncol(binned_mat),
                       dimnames = dimnames(binned_mat))

cell_order <- cell_clust$order
ordered_cells <- colnames(binned_mat)[cell_order]
ordered_clusters <- cell_clusters[cell_order]

for (cl in unique(ordered_clusters)) {
  cluster_cells <- ordered_cells[ordered_clusters == cl]
  cluster_mat <- binned_mat[, cluster_cells, drop = FALSE]
  n_cells <- ncol(cluster_mat)
  
  # Apply rolling mean row-wise (i.e., across cells per bin)
  smoothed_cluster_mat <- matrix(NA, nrow = nrow(cluster_mat), ncol = n_cells)
  
  for (i in 1:nrow(cluster_mat)) {
    bin_row <- cluster_mat[i, ]
    smoothed <- numeric(n_cells)
    
    for (j in seq_len(n_cells)) {
      left <- max(1, j - half_window)
      right <- min(n_cells, j + half_window)
      smoothed[j] <- mean(bin_row[left:right], na.rm = TRUE)
    }
    
    smoothed_cluster_mat[i, ] <- smoothed
  }
  
  smoothed_mat[, cluster_cells] <- smoothed_cluster_mat
}

# pdf(paste0("CNV_results/", unique(meta$ident), "_CNV_profile_smoothed.pdf"), width = 15, height = 8)
# draw(
  Heatmap(
    smoothed_mat,
    name = "CNV",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_order = cell_clust$order,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = ha,
    column_title = "CNV Heatmap (Binned by Genomic Position)",
    row_title_rot = 0,
    use_raster = FALSE
  )
#)
# dev.off()


binned_mat_colored <- matrix(0, nrow = nrow(smoothed_mat), ncol = ncol(smoothed_mat),
                             dimnames = dimnames(smoothed_mat))
binned_mat_colored[smoothed_mat > 0.0] <- 1
binned_mat_colored[smoothed_mat < -0.0] <- -1
col_fun <- c("-1" = "blue", "0" = "white", "1" = "red")


# Distance matrix and hierarchical clustering on transposed matrix (cells)
cell_dist <- dist(t(binned_mat_colored))  # Cells are columns
cell_clust <- hclust(cell_dist, method = "ward.D2")

# Cut tree into desired number of clusters (adjust k as needed)
k <- 4
cell_clusters <- cutree(cell_clust, k = k)

# Add cluster annotation to metadata
meta_filtered$cluster <- factor(cell_clusters[colnames(binned_mat_colored)])

# Assign colors to clusters
cluster_ids <- sort(unique(meta_filtered$cluster))

coord <- cnaScatterPlot(outs)
meta_filtered$cna_signal <- coord$cna.signal
meta_filtered$cna_cor <- coord$cna.cor

valid_cluster_ids <- sort(unique(na.omit(meta_filtered$cluster)))

n_colors <- max(3, min(length(valid_cluster_ids), 8))
cluster_colors <- setNames(
  RColorBrewer::brewer.pal(n_colors, "Set2")[seq_along(valid_cluster_ids)],
  valid_cluster_ids
)
valid_idents <- sort(unique(na.omit(meta_filtered$ident)))

ident_colors <- ident_colors[valid_idents]
ident_colors <- setNames(ident_colors, valid_idents)

aggregate(cbind(cna_signal, cna_cor) ~ cluster, data = meta_filtered, FUN = mean)

cluster_mean_signal <- aggregate(cna_signal ~ cluster, data = meta_filtered, FUN = mean)
cluster_mean_cor <- aggregate(cna_cor ~ cluster, data = meta_filtered, FUN = mean)

# Create named vectors for easy mapping
mean_signal_map <- setNames(cluster_mean_signal$cna_signal, cluster_mean_signal$cluster)
mean_cor_map <- setNames(cluster_mean_cor$cna_cor, cluster_mean_cor$cluster)

# Map values to all cells in meta_filtered
meta_filtered$mean_signal_per_cluster <- mean_signal_map[as.character(meta_filtered$cluster)]
meta_filtered$mean_cor_per_cluster <- mean_cor_map[as.character(meta_filtered$cluster)]

library(circlize)

max_signal <- max(meta_filtered$mean_signal_per_cluster, na.rm = TRUE)
max_cor <- max(meta_filtered$mean_cor_per_cluster, na.rm = TRUE)

col_list <- list(
  orig.ident = ident_colors,
  cluster = cluster_colors,
  mean_cna_signal = colorRamp2(c(0, max_signal), c("white", "red")),
  mean_cna_cor = colorRamp2(c(0, max_cor), c("white", "blue"))
)

# Categorize cna_cor values
meta_filtered$cna_cor_cat <- cut(
  meta_filtered$cna_cor,
  breaks = c(-Inf, 0.15, 0.25, Inf),
  labels = c("<0.15", "0.15–0.25", ">0.25"),
  right = FALSE
)

# Define colors for cna_cor_cat
cna_cor_cat_colors <- setNames(
  c("#D3D3D3", "#87CEFA", "#4169E1"),
  c("<0.15", "0.15–0.25", ">0.25")
)

# Add this to the annotation color list
col_list$cna_cor_cat <- cna_cor_cat_colors

# Update annotation block
ha <- HeatmapAnnotation(
  orig.ident = meta_filtered$ident,
  cluster = meta_filtered$cluster,
  mean_cna_signal = meta_filtered$mean_signal_per_cluster,
  mean_cna_cor = meta_filtered$mean_cor_per_cluster,
  cna_cor_cat = meta_filtered$cna_cor_cat,
  col = col_list,
  annotation_name_side = "left",
  annotation_legend_param = list(
    cluster = list(title = "Cluster"),
    orig.ident = list(title = "Sample"),
    mean_cna_signal = list(title = "Mean CNA Signal"),
    mean_cna_cor = list(title = "Mean CNA Correlation"),
    cna_cor_cat = list(title = "CNA Correlation Category")
  ),
  annotation_height = unit(c(5, 5, 4, 4, 4), "mm"),
  simple_anno_size = unit(2, "mm")
)

# pdf(paste0("CNV_results/", unique(meta$ident), "_CNV_profile_binary.pdf"), width = 15, height = 8)
Heatmap(
  binned_mat_colored,
  name = "CNV",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = cell_clust$order,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  column_title = "CNV Heatmap (Binned, Thresholded)",
  row_title_rot = 0,
  use_raster = FALSE,
  col = col_fun
)
# dev.off()



low_clust <- cluster_mean_cor$cluster[cluster_mean_cor$cna_cor < 0.1]
selected_cells <- colnames(binned_mat_colored)[
  (meta_filtered$cna_cor < 0.25 | meta_filtered$cluster %in% low_clust) & meta_filtered$ident == "patient_H"
]
mat_subset <- binned_mat_colored[, !(colnames(binned_mat_colored) %in% selected_cells), drop = FALSE]
meta_subset <- meta_filtered[!(colnames(binned_mat_colored) %in% selected_cells), , drop = FALSE]

ha_subset <- HeatmapAnnotation(
  orig.ident = meta_subset$ident,
  cluster = meta_subset$cluster,
  mean_cna_signal = meta_subset$mean_signal_per_cluster,
  mean_cna_cor = meta_subset$mean_cor_per_cluster,
  cna_cor_cat = meta_subset$cna_cor_cat,
  col = col_list,
  annotation_name_side = "left",
  annotation_legend_param = list(
    cluster = list(title = "Cluster"),
    orig.ident = list(title = "Sample"),
    mean_cna_signal = list(title = "Mean CNA Signal"),
    mean_cna_cor = list(title = "Mean CNA Correlation"),
    cna_cor_cat = list(title = "CNA Correlation Category")
  ),
  annotation_height = unit(c(5, 4), "mm"),
  simple_anno_size = unit(2, "mm")
)

ordered_cells <- colnames(binned_mat_colored)[cell_clust$order]
ordered_subset <- ordered_cells[!(ordered_cells %in% selected_cells)]

pdf(paste0("CNV_results/", unique(meta$ident), "_CNV_profile_binary_filtered.pdf"), width = 15, height = 8)
Heatmap(
  mat_subset,
  name = "CNV",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = ordered_subset,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_subset,
  column_title = "",
  row_title_rot = 0,
  use_raster = FALSE,
  col = col_fun
)
dev.off()


orig_mat <- binned_mat[, colnames(mat_subset)]
norm_mat <- sweep(orig_mat, 2, colSums(abs(orig_mat)), FUN = "/")
norm_mat <- sign(norm_mat) * (abs(norm_mat) ^ 1.5) * 100

Heatmap(
  norm_mat,
  name = "CNV",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = ordered_subset, 
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_subset,
  column_title = "",
  row_title_rot = 0,
  use_raster = FALSE
)
