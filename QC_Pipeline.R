###############Loding required packages########################

library("Seurat")
library("dplyr")
library("Seurat")
library("patchwork")
library("ggplot2")
library("foreach")
library("doParallel")
library("circlize")
library("gridExtra")
library("grid")
library("tidyr")
library("tibble")
library("purrr")
library("readxl")
library('parallel')
library('future')

setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

######################Setting parameters########################

batch <- "_PDOs"
names_tmdata <- list.files(path = "/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all", full.names = FALSE, recursive = FALSE)
names_tmdata <- sub(".csv", "", names_tmdata[grepl("*_PDO.csv", names_tmdata)])
n_clusters = 8

qc_rules <- data.frame(
  pattern = c(
    "_Untreated",
    "_Treated",
    "_PDO"
  ),
  mito   = c(15, 15, 15),
  nGenes_min = c(500, 500, 500),
  nGenes_max = c(7000, 7000, 13000), 
  hk     = c(3, 3, 3),
  stringsAsFactors = FALSE
)

################################################################

initialise <- function(names) {
  
  tmdata_list <- list()
  for (name in names) {
    filename <- paste0("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all/", name, ".csv")
    tmdata <- data.table::fread(filename)
    Genes <- tmdata[[1]]
    counts <- as.matrix(tmdata[, -1])
    rownames(counts) <- Genes
    colnames(counts) <- make.unique(colnames(counts))
    rownames(counts) <- make.unique(rownames(counts))
    tmdata_list[[name]] <- CreateSeuratObject(counts = counts)
    tmdata_list[[name]]$orig.ident <- rep(name, dim(tmdata_list[[name]])[2])
    print(paste0("finished reading ", name))
  }
  
  return(tmdata_list)
}

inspect <- function(tmdata_list) {
  
  x_features_plot <- list()
  x_count_plot <- list()
  x_mito_plot <- list()
  
  for (name in names(tmdata_list)) {
    
    tmdata_list[[name]][["percent.mt"]] <- PercentageFeatureSet(tmdata_list[[name]], pattern = "^MT-")
    
    mean_nFeature <- mean(tmdata_list[[name]]$nFeature_RNA, na.rm = TRUE)
    median_nFeature <- median(tmdata_list[[name]]$nFeature_RNA, na.rm = TRUE)
    
    mean_nCount <- mean(tmdata_list[[name]]$nCount_RNA, na.rm = TRUE)
    median_nCount <- median(tmdata_list[[name]]$nCount_RNA, na.rm = TRUE)
    
    mean_percent_mt <- mean(tmdata_list[[name]]$percent.mt, na.rm = TRUE)
    median_percent_mt <- median(tmdata_list[[name]]$percent.mt, na.rm = TRUE)
    
    base_theme <- theme(
      text = element_text(size = 8),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 6),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 6),
      legend.position = "none"
    )
    
    # nFeature_RNA plot
    x_features_plot[[name]] <- VlnPlot(tmdata_list[[name]], features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") +
      base_theme +
      geom_hline(yintercept = mean_nFeature, linetype = "dashed", color = "blue", size = 0.5) +
      geom_hline(yintercept = median_nFeature, linetype = "solid", color = "red", size = 0.5) +
      annotate("text", x = 1.5, y = mean_nFeature, label = paste("Mean:", round(mean_nFeature, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1.5, y = median_nFeature, label = paste("Median:", round(median_nFeature, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(tmdata_list[[name]])),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    
    # nCount_RNA plot
    x_count_plot[[name]] <- VlnPlot(tmdata_list[[name]], features = "nCount_RNA", pt.size = 0, group.by = "orig.ident") +
      base_theme +
      geom_hline(yintercept = mean_nCount, linetype = "dashed", color = "blue", size = 0.5) +
      geom_hline(yintercept = median_nCount, linetype = "solid", color = "red", size = 0.5) +
      annotate("text", x = 1.5, y = mean_nCount, label = paste("Mean:", round(mean_nCount, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1.5, y = median_nCount, label = paste("Median:", round(median_nCount, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(tmdata_list[[name]])),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black")
    
    # percent.mt plot
    x_mito_plot[[name]] <- VlnPlot(tmdata_list[[name]], features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
      base_theme +
      geom_hline(yintercept = mean_percent_mt, linetype = "dashed", color = "blue", size = 0.5) +
      geom_hline(yintercept = median_percent_mt, linetype = "solid", color = "red", size = 0.5) +
      annotate("text", x = 1.5, y = mean_percent_mt, label = paste("Mean:", round(mean_percent_mt, 1)),
               hjust = 0.5, vjust = -1, size = 3, color = "blue") +
      annotate("text", x = 1.5, y = median_percent_mt, label = paste("Median:", round(median_percent_mt, 1)),
               hjust = 0.5, vjust = 1.5, size = 3, color = "red") +
      annotate("text", x = Inf, y = Inf, label = paste("NCells:", ncol(tmdata_list[[name]])),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black")
  }
  
  plot_chunks <- function(plot_list) {
    split(plot_list, ceiling(seq_along(plot_list) / 6))
  }
  
  get_grid_dims <- function(n) {
    if (n == 1) return(c(1, 1))
    if (n == 2) return(c(1, 2))
    if (n == 3) return(c(2, 2))
    if (n == 4) return(c(2, 2))
    return(c(2, 3)) # Default for larger sets
  }
  
  pdf(paste0("Inspections", batch, ".pdf"), width = 8, height = 11)
  
  for (plot_list in list(x_features_plot, x_count_plot, x_mito_plot)) {
    for (chunk in plot_chunks(plot_list)) {
      dims <- get_grid_dims(length(chunk))
      grid.arrange(grobs = chunk, ncol = dims[1], nrow = dims[2])
    }
  }
  
  dev.off()
}

normalise <- function(tmdata_list) {
  
  for (name in names(tmdata_list)) {
    CPM <- apply(tmdata_list[[name]]@assays$RNA$counts, 2, function(x) (x / sum(x)) * 1e6)
    CPM <- as(CPM, "dgCMatrix")
    tmdata_list[[name]]@assays$RNA$CPM <- CPM
    expr <- log2((CPM / 10) + 1)
    expr <- as(expr, "CsparseMatrix")
    tmdata_list[[name]]@assays$RNA$data <- expr
    print(paste0("finished normalising ", name))
  }
  
  return(tmdata_list)
}

doublets_filtering <- function(tmdata, sampleid) {
  
  tmdata <- NormalizeData(tmdata)
  tmdata <- FindVariableFeatures(tmdata, selection.method = "vst", nfeatures = 2000)
  tmdata <- ScaleData(tmdata)
  tmdata <- RunPCA(tmdata)
  tmdata <- FindNeighbors(tmdata, dims = 1:50)
  tmdata <- FindClusters(tmdata, resolution = 0.5)
  tmdata <- RunUMAP(tmdata, dims = 1:50)
  
  options(future.globals.maxSize = 2 * 1024^3)
  sweep_tmdata <- paramSweep(tmdata, PCs = 1:50, sct = T)
  sweep_summary <- summarizeSweep(sweep_tmdata)
  bcmvn <- find.pK(sweep_summary)
  pk_tmdata <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  
  ncells <- ncol(tmdata)
  if (ncells <= 1000) {
    DR <- 0.008
  } else if (ncells <= 5000) {
    DR <- 0.04
  } else if (ncells <= 10000) {
    DR <- 0.08
  } else if (ncells <= 20000) {
    DR <- 0.16
  } else if (ncells <= 30000) {
    DR <- 0.24
  } else {
    DR <- 0.24
  }
  homotypic <- modelHomotypic(tmdata$seurat_clusters)
  nExp <- round(DR * ncol(tmdata))
  nExp_adj <- round(nExp * (1 - homotypic))
  
  tmdata <- doubletFinder(tmdata, PCs = 1:50, pN = 0.25, pK = pk_tmdata,
                          nExp = nExp_adj, reuse.pANN = F, sct = F)
  
  SorD <- grep("DF", colnames(tmdata@meta.data), value = T)
  p1 <- DimPlot(tmdata, reduction = "umap")
  p2 <- DimPlot(tmdata, reduction = "umap", group.by = SorD)
  combined_plot <- p1 + p2
  
  singlet <- colnames(tmdata)[which(tmdata@meta.data[[SorD]] == "Singlet")]
  tmdata <- subset(tmdata, cells = singlet)
  
  return(list(tmdata = tmdata, plot = combined_plot))
}

doublets_parallel <- function(tmdata_list) {
  
  cl <- makeCluster(n_clusters)
  registerDoParallel(cl)
  
  clusterExport(cl, c("doublets_filtering"))
  
  tmdata_updated <- foreach(name = names(tmdata_list),
                            .packages = c("Seurat", "DoubletFinder", "ggplot2", "future"),
                            .errorhandling = 'pass') %dopar% {
                              tryCatch({
                                doublets_filtering(
                                  tmdata = tmdata_list[[name]],
                                  sampleid = name
                                )
                              }, error = function(e) {
                                message("Error in sample ", name, ": ", conditionMessage(e))
                                list(tmdata = tmdata_list[[name]], plot = NULL, error = conditionMessage(e))
                              })
                            }
  stopCluster(cl)
  names(tmdata_updated) <- names(tmdata_list)
  
  pdf(paste0("doublets_filtering", batch, ".pdf"), width = 12, height = 8)
  for (name in names(tmdata_updated)) {
    if (!is.null(tmdata_updated[[name]]$plot)) {
      print(tmdata_updated[[name]]$plot + plot_annotation(title = name))
    }
  }
  dev.off()
  
  tmdata_output <- lapply(tmdata_updated, function(x) x$tmdata)
  names(tmdata_output) <- names(tmdata_updated)

  error_messages <- sapply(names(tmdata_updated), function(nm) tmdata_updated[[nm]]$error)
  error_messages <- error_messages[!sapply(error_messages, is.null)]
  if (length(error_messages) > 0) {
    cat("\nSummary of errors:\n")
    print(error_messages)
  }
  
  return(tmdata_output)
}

cells_filtering <- function(tmdata_list, rules = qc_rules) {
  
  plot <- list()
  g_filter <- vector()
  hk_filter <- vector()
  for (name in names(tmdata_list)) {
    ##############
    match_idx <- which(sapply(qc_rules$pattern, function(p) grepl(p, name)))
    matched_row <- qc_rules[match_idx[1], ]
    ngenes_min <- matched_row$nGenes_min
    ngenes_max <- matched_row$nGenes_max
    hkmean <- matched_row$hk
    ##############
    expr <- as.matrix(tmdata_list[[name]]@assays$RNA$data)
    n_genes <- colSums(expr > 0)
    hk_list <- c("ACTB", "GAPDH", "RPS11", "RPS13", "RPS14", "RPS15", "RPS16", "RPS18",
                 "RPS19", "RPS20", "RPL10", "RPL13", "RPL15", "RPL18")
    
    hk_list <- hk_list[hk_list %in% rownames(expr)]
    hk_expression <- expr[hk_list, , drop = FALSE]
    hk_mean <- colMeans(hk_expression)
    sl_cells_g <- n_genes >= ngenes_min & n_genes <= ngenes_max
    g_filter[[name]] <- sum(sl_cells_g)
    sl_cells_hk <- hk_mean >= hkmean & sl_cells_g
    hk_filter[[name]] <- sum(sl_cells_hk)
    if (sum(sl_cells_hk) != 0) {
      tmdata_list[[name]] <- subset(tmdata_list[[name]], cells = names(sl_cells_hk)[sl_cells_hk])
    } else {
      tmdata_list[[name]] <- NULL
    }
    
    
    plot_data <- data.frame(hk_mean = hk_mean, n_genes = n_genes, valid = sl_cells_hk)
    
    p <- ggplot(plot_data, aes(x = n_genes, y = hk_mean, color = valid)) +
      geom_point() +
      scale_x_continuous(trans = "log10", labels = scales::comma) +
      scale_y_continuous(trans = "log10", labels = scales::comma) +
      scale_color_manual(values = c("lightgrey", "black")) +
      labs(x = "Number of Genes", y = "HK Mean", color = "Valid Cells") +
      theme_minimal() +
      theme(legend.position = "none", plot.title = element_text(size = 8)) +
      annotate("text", x = Inf, y = Inf, label = paste0("NCells passed: ", sum(sl_cells_hk)),
               hjust = 1.1, vjust = 1.1, size = 3, color = "black") +
      ggtitle(name) +
      geom_vline(xintercept = ngenes_min, linetype = "dashed", color = "red") + 
      geom_vline(xintercept = ngenes_max, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = hkmean, linetype = "dashed", color = "red") + 
      annotate("text", x = ngenes_min, y = min(plot_data$hk_mean), label = paste0("Number of genes > ", ngenes_min),
               hjust = -0.1, vjust = 0, size = 3, color = "red") + 
      annotate("text", x = ngenes_max, y = min(plot_data$hk_mean), label = paste0("Number of genes < ", ngenes_max),
               hjust = 1.1, vjust = 0, size = 3, color = "red") + 
      annotate("text", x = min(plot_data$n_genes), y = hkmean, label = paste0("HK Mean > ", hkmean),
               hjust = 0, vjust = -0.5, size = 3, color = "red")
    
    plot[[name]] <- p
    print(paste0("finished cell filtering for  ", name))
  }
  
  split_plots <- function(plot_list) {
    split(plot_list, ceiling(seq_along(plot_list) / 6))
  }
  get_layout_dims <- function(n) {
    if (n == 1) return(c(1, 1))
    if (n == 2) return(c(1, 2))
    if (n == 3 || n == 4) return(c(2, 2))
    return(c(2, 3))  # Default layout for 5-6 plots
  }
  
  pdf(paste0("cells_filtering", batch, ".pdf"), width = 8, height = 11)
  
  for (chunk in split_plots(plot)) {
    dims <- get_layout_dims(length(chunk))
    grid.arrange(grobs = chunk, ncol = dims[1], nrow = dims[2])  # ncol and nrow reversed
  }
  dev.off()
  
  return(list(tmdata_list, g_filter, hk_filter))
}

###################################################################

raw <- initialise(names_tmdata)
x_filter <- sapply(raw, ncol)
filtered <- raw
rm(raw)
inspect(filtered)
filtered <- doublets_parallel(filtered)
print("finished finding doublets")
db_filter <- sapply(filtered, ncol)
for (name in names(filtered)) {
  filtered[[name]][["percent.mt"]] <- PercentageFeatureSet(filtered[[name]], pattern = "^MT-")
  match_idx <- which(sapply(qc_rules$pattern, function(p) grepl(p, name)))
  matched_row <- qc_rules[match_idx[1], ]
  max_mt <- matched_row$mito
  filtered[[name]] <- subset(filtered[[name]], percent.mt < max_mt)
}
mt_filter <- sapply(filtered, ncol)
filtered <- normalise(filtered)
cells_ft_outs <- cells_filtering(filtered, rules = qc_rules)
filtered <- cells_ft_outs[[1]]
g_filter <- cells_ft_outs[[2]]
hk_filter <- cells_ft_outs[[3]]

sm_table <- cbind(x_filter, db_filter, mt_filter, g_filter, hk_filter)
colnames(sm_table) <- c("raw", "doublets\nremoved", paste0("mito_DNA\npercentage < 15"),
                        paste0("number of\ngenes filter"), paste0("housekeeping\nexpression > 3"))

write.csv(sm_table, paste0("filtering_summary", batch, ".csv"))

saveRDS(filtered, paste0("PDOs_list", batch, ".rds"))

out_dir <- "by_samples"; if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
mclapply(names(filtered), function(nm) {
  sample_dir <- file.path(out_dir, nm)
  if (!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
  saveRDS(filtered[[nm]], file.path(sample_dir, paste0(nm, ".rds")), compress = FALSE)
  NULL
}, mc.cores = n_clusters, mc.preschedule = FALSE)

###############################################################################

all_genes <- Reduce(union, lapply(filtered, function(obj) {
  rownames(GetAssayData(obj, layer = "counts"))
}))

pad_matrix <- function(mat, all_genes) {
  missing_genes <- setdiff(all_genes, rownames(mat))
  if (length(missing_genes) > 0) {
    zero_mat <- Matrix::Matrix(0, nrow = length(missing_genes), ncol = ncol(mat),
                               dimnames = list(missing_genes, colnames(mat)))
    mat <- rbind(mat, zero_mat)
  }
  # Ensure the same row order
  mat <- mat[all_genes, , drop = FALSE]
  return(mat)
}

counts_list <- lapply(names(filtered), function(id) {
  mat <- GetAssayData(filtered[[id]], layer = "counts")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})

cpm_list <- lapply(names(filtered), function(id) {
  mat <- GetAssayData(filtered[[id]], layer = "CPM")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})

lognorm_list <- lapply(names(filtered), function(id) {
  mat <- GetAssayData(filtered[[id]], layer = "data")
  colnames(mat) <- paste(id, colnames(mat), sep = "_")
  pad_matrix(mat, all_genes)
})

meta_list <- lapply(names(filtered), function(id) {
  meta <- filtered[[id]]@meta.data[, c(names(filtered[[id]]@meta.data)[1:4])]
  rownames(meta) <- paste(id, rownames(meta), sep = "_")
  return(meta)
})

combined_counts <- do.call(cbind, counts_list)
combined_cpm <- do.call(cbind, cpm_list)
combined_lognorm <- do.call(cbind, lognorm_list)
combined_meta <- do.call(rbind, meta_list)

merged_obj <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta)

merged_obj@assays$RNA$CPM <- combined_cpm
merged_obj@assays$RNA$data <- combined_lognorm
rm(filtered)

merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
merged_obj <- FindNeighbors(merged_obj, dims = 1:50)

merged_obj <- FindClusters(merged_obj, resolution = 0.8, algorithm = 1)
merged_obj$leiden_clusters <- Idents(merged_obj)

merged_obj <- RunUMAP(merged_obj, dims = 1:50)

library(readxl)
data <- read_excel("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/SP_Nicola work_amended_michael_Keito-190825.xlsx", sheet = 1)
data <- as.matrix(data)
rownames(data) <- data[, 1]
data <- t(data[ , -1])
rownames(data) <- paste0("SUR", rownames(data))
data <- data.frame(SUR = rownames(data), data, row.names = NULL)

merged_obj$SUR <- sub("^(SUR[0-9]+).*", "\\1", merged_obj$orig.ident)
merged_obj$Batch <- sub("^[^_]+_", "", merged_obj$orig.ident)
cell_names <- colnames(merged_obj)
merged_obj@meta.data <- merged_obj@meta.data %>%
  left_join(data, by = "SUR")
rownames(merged_obj@meta.data) <- cell_names

saveRDS(merged_obj, paste0("PDOs_merged.rds"))
saveRDS(merged_obj@meta.data, "PDOs_all_meta.rds")

p1 <- DimPlot(merged_obj, group.by = "leiden_clusters", label = TRUE) + ggtitle("Louvain Clustering")
p2 <- DimPlot(merged_obj, group.by = "orig.ident", label = FALSE) + ggtitle("orig.ident")
p3 <- DimPlot(merged_obj, group.by = "Batch", label = FALSE) + ggtitle("Batch")
p4 <- DimPlot(merged_obj, group.by = "Clinical.response.at.OG.MDT..responder.non.responder", label = FALSE) + ggtitle("Response")

combined_plot <- (p1 | p2) / (p3 | p4)

ggsave("PDOs_merged.png", plot = combined_plot, width = 12, height = 6, dpi = 300)
