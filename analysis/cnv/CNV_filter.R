library(MASS)      # fitdistr
library(dplyr)
library(tidyr)
library(matrixStats)

gene_order <- read.table("/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt", header = FALSE)
colnames(gene_order) <- c("gene_id", "chromosome", "start", "end")
common_genes <- intersect(rownames(outs), gene_order$gene_id)
outs <- outs[common_genes, ]
gene_order <- gene_order %>%
  filter(gene_id %in% rownames(outs))

gene_order <- gene_order[match(rownames(outs), gene_order$gene_id), ]

genes_by_chr <- split(gene_order$gene_id, gene_order$chromosome)


########################################################################

cna_signal <- vapply(
  genes_by_chr,
  function(g) {
    if (length(g) == 1) {
      outs[g, ]                       # keep as numeric vector
    } else {
      colMeans(outs[g, , drop = FALSE])
    }
  },
  numeric(ncol(outs))
)

cna_signal <- t(cna_signal)
chromosomes <- rownames(cna_signal)

ref_cells <- ref
ref_cells <- unlist(ref_cells, use.names = FALSE)

fit_params <- lapply(chromosomes, function(chr) {
  vals <- cna_signal[chr, ref_cells]
  fit  <- fitdistr(vals, "normal")
  c(mean = fit$estimate["mean"], sd = fit$estimate["sd"])
})
names(fit_params) <- chromosomes



zmat   <- matrix(NA_real_, nrow = nrow(cna_signal), ncol = ncol(cna_signal),
                 dimnames = dimnames(cna_signal))
p_mat  <- zmat
p_adj  <- zmat

for (chr in chromosomes) {
  mu <- fit_params[[chr]]["mean.mean"]
  sd <- fit_params[[chr]]["sd.sd"]
  z  <- (cna_signal[chr, ] - mu) / sd
  p  <- 2 * pnorm(-abs(z))                     # two-sided
  
  zmat[chr, ]  <- z
  p_mat[chr, ] <- p
  p_adj[chr, ] <- p.adjust(p, method = "BH")
}



call_matrix <- matrix(NA_character_, nrow = nrow(cna_signal), 
                      ncol = ncol(cna_signal), dimnames = dimnames(cna_signal))

is_event <- p_adj < 0.05
is_susp  <- !is_event & (p_mat < 0.05)

call_matrix[is_event & (cna_signal > 0)] <- "AMP"
call_matrix[is_event & (cna_signal < 0)] <- "DEL"

call_matrix[is_susp  & (cna_signal > 0)] <- "susp_AMP"
call_matrix[is_susp  & (cna_signal < 0)] <- "susp_DEL"



cna_data <- read.delim("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/matched/CNA_Genes.txt", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cna_data <- cna_data[cna_data$Profiled.Samples > 30, ]

panel_amp_chr <- unique(gene_order$chr[gene_order$gene %in% 
                                         cna_data$Gene[cna_data$CNA == "AMP"]])
panel_del_chr <- unique(gene_order$chr[gene_order$gene %in% 
                                         cna_data$Gene[cna_data$CNA == "HOMDEL"]])

for (chr in chromosomes) {
  is_amp <- call_matrix[chr, ] == "AMP" & chr %in% panel_amp_chr
  is_del <- call_matrix[chr, ] == "DEL" & chr %in% panel_del_chr
  
  call_matrix[chr, is_amp] <- "real_AMP"
  call_matrix[chr, is_del] <- "real_DEL"
}


celltype_init <- meta[[ct_col]][match(colnames(outs), rownames(meta))]
names(celltype_init) <- colnames(outs)

cell_status <- rep(NA_character_, length(celltype_init))
names(cell_status) <- names(celltype_init)

has_real <- apply(call_matrix, 2, function(col) any(grepl("^real_", col)))
has_suspicious <- apply(call_matrix, 2, function(col) any(grepl("^susp_", col)))

for (id in names(cell_status)) {
  is_ref <- id %in% ref_cells
  
  if (is_ref) {                                 # originally non-malignant
    if (has_real[id] | has_suspicious[id]) {
      cell_status[id] <- "unresolved"
    } else {
      cell_status[id] <- "nonmalignant"
    }
    
  } else {                                      # presumed malignant
    if (has_real[id]) {
      cell_status[id] <- "malignant"
    } else if (has_suspicious[id]) {
      cell_status[id] <- "unresolved"
    } else {
      cell_status[id] <- "nonmalignant"
    }
  }
}

event_long <- as.data.frame(call_matrix) %>% 
  tibble::rownames_to_column("chromosome") %>% 
  pivot_longer(-chromosome, names_to = "cell_id", values_to = "call")

cell_summary <- data.frame(
  cell_id   = names(cell_status),
  ident = meta$ident, 
  classification = cell_status,
  celltypes = meta[[ct_col]], 
  stringsAsFactors = FALSE
)


cell_status[cell_summary$classification == "malignant"]

########################################################################


cna_signal <- outs  # rows = genes, cols = cells
rownames(cna_signal) <- unlist(genes_by_chr, use.names = FALSE)

# Sort gene_order and cna_signal together
gene_order <- gene_order[match(rownames(cna_signal), gene_order$gene_id), ]
stopifnot(all.equal(rownames(cna_signal), gene_order$gene_id))

# Initialize result matrices
zmat   <- matrix(NA_real_, nrow = nrow(cna_signal), ncol = ncol(cna_signal), dimnames = dimnames(cna_signal))
p_mat  <- zmat
p_adj  <- zmat

ref_cells <- unlist(ref, use.names = FALSE)
window_size <- 50

# Apply sliding window per chromosome
for (chr in unique(gene_order$chromosome)) {
  gene_ids_chr <- gene_order$gene_id[gene_order$chromosome == chr]
  chr_idx <- which(rownames(cna_signal) %in% gene_ids_chr)
  num_genes_chr <- length(chr_idx)
  
  for (i in seq_along(chr_idx)) {
    global_idx <- chr_idx[i]
    rel_idx <- i
    
    start_idx <- max(1, rel_idx - floor(window_size / 2))
    end_idx   <- min(num_genes_chr, rel_idx + floor(window_size / 2))
    window_genes <- chr_idx[start_idx:end_idx]
    
    ref_vals <- cna_signal[window_genes, ref_cells, drop = FALSE]
    flat_vals <- as.numeric(ref_vals)
    
    fit <- tryCatch(
      fitdistr(flat_vals, "normal"),
      error = function(e) list(estimate = c(mean = mean(flat_vals), sd = sd(flat_vals)))
    )
    
    mu <- fit$estimate["mean"]
    sd <- fit$estimate["sd"]
    gene_vals <- cna_signal[global_idx, ]
    
    z <- (gene_vals - mu) / sd
    p <- 2 * pnorm(-abs(z))  # two-sided
    
    zmat[global_idx, ] <- z
    p_mat[global_idx, ] <- p
    p_adj[global_idx, ] <- p.adjust(p, method = "BH")
  }
}

# Initial CNA calls
call_matrix <- matrix(NA_character_, nrow = nrow(cna_signal), ncol = ncol(cna_signal), dimnames = dimnames(cna_signal))
is_event <- p_adj < 0.05
is_susp  <- !is_event & (p_mat < 0.05)

call_matrix[is_event & (cna_signal > 0)] <- "AMP"
call_matrix[is_event & (cna_signal < 0)] <- "DEL"
call_matrix[is_susp  & (cna_signal > 0)] <- "susp_AMP"
call_matrix[is_susp  & (cna_signal < 0)] <- "susp_DEL"

# Load panel CNA data
cna_data <- read.delim("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/matched/CNA_Genes.txt", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cna_data <- cna_data[cna_data$Profiled.Samples > 36, ]

panel_amp_genes <- intersect(cna_data$Gene[cna_data$CNA == "AMP"], rownames(cna_signal))
panel_del_genes <- intersect(cna_data$Gene[cna_data$CNA == "HOMDEL"], rownames(cna_signal))

# Function to get gene index window (within the same chromosome)
get_window_indices <- function(gene, gene_order_chr_idx) {
  gene_pos <- which(gene_order_chr_idx == gene)
  start <- max(1, gene_pos - floor(window_size / 2))
  end <- min(length(gene_order_chr_idx), gene_pos + floor(window_size / 2))
  which(rownames(cna_signal) %in% gene_order_chr_idx[start:end])
}

# Reassign to real_AMP / real_DEL if matched to panel genes
for (chr in unique(gene_order$chromosome)) {
  gene_order_chr_idx <- gene_order$gene_id[gene_order$chromosome == chr]
  
  amp_genes_chr <- intersect(panel_amp_genes, gene_order_chr_idx)
  del_genes_chr <- intersect(panel_del_genes, gene_order_chr_idx)
  
  for (gene in amp_genes_chr) {
    idx <- get_window_indices(gene, gene_order_chr_idx)
    call_matrix[idx, ][call_matrix[idx, ] == "AMP"] <- "real_AMP"
  }
  for (gene in del_genes_chr) {
    idx <- get_window_indices(gene, gene_order_chr_idx)
    call_matrix[idx, ][call_matrix[idx, ] == "DEL"] <- "real_DEL"
  }
}


# Count real AMP or DEL per cell (column-wise)
real_event_counts <- colSums((call_matrix == "real_AMP" | call_matrix == "real_DEL") & !is.na(call_matrix))

# Optionally, name the vector
names(real_event_counts) <- colnames(call_matrix)

celltype_init <- meta[[ct_col]][match(colnames(outs), rownames(meta))]
names(celltype_init) <- colnames(outs)

cell_status <- rep(NA_character_, length(celltype_init))
names(cell_status) <- names(celltype_init)

has_real <- apply(call_matrix, 2, function(col) any(grepl("^real_", col)))
has_suspicious <- apply(call_matrix, 2, function(col) any(grepl("^susp_", col)))

for (id in names(cell_status)) {
  is_ref <- id %in% ref_cells
  
  if (is_ref) {                                 # originally non-malignant
    if (has_real[id] | has_suspicious[id]) {
      cell_status[id] <- "unresolved"
    } else {
      cell_status[id] <- "nonmalignant"
    }
    
  } else {                                      # presumed malignant
    if (has_real[id]) {
      cell_status[id] <- "malignant"
    } else if (has_suspicious[id]) {
      cell_status[id] <- "unresolved"
    } else {
      cell_status[id] <- "nonmalignant"
    }
  }
}

event_long <- as.data.frame(call_matrix) %>% 
  tibble::rownames_to_column("chromosome") %>% 
  pivot_longer(-chromosome, names_to = "cell_id", values_to = "call")

cell_summary <- data.frame(
  cell_id   = names(cell_status),
  ident = meta$ident, 
  classification = cell_status,
  celltypes = meta[[ct_col]], 
  stringsAsFactors = FALSE
)

cell_status[cell_summary$classification == "malignant"]

table(cell_summary[cell_summary$ident == "patient_L", ]$classification)
table(cell_summary[cell_summary$ident == "patient_L" & cell_summary$celltypes == "Ductal", ]$classification)
table(cell_summary[cell_summary$ident == "sn_reference", ]$classification)

table(cell_summary[cell_summary$ident == "SUR791", ]$classification)
table(cell_summary[cell_summary$ident == "SUR791" & cell_summary$celltypes == "Ductal", ]$classification)
table(cell_summary[cell_summary$ident == "sc_reference", ]$classification)

table(cell_summary[cell_summary$ident == "D02_Normal", ]$classification)
table(cell_summary[cell_summary$ident == "D02_Normal" & cell_summary$celltypes == "Epithelial", ]$classification)

table(cell_summary[cell_summary$ident == "D03_Cancer", ]$classification)
table(cell_summary[cell_summary$ident == "D03_Cancer" & cell_summary$celltypes == "Epithelial", ]$classification)
table(cell_summary[cell_summary$ident == "D04_Cancer", ]$classification)
table(cell_summary[cell_summary$ident == "D04_Cancer" & cell_summary$celltypes == "Epithelial", ]$classification)
table(cell_summary[cell_summary$ident == "D05_Cancer", ]$classification)
table(cell_summary[cell_summary$ident == "D05_Cancer" & cell_summary$celltypes == "Epithelial", ]$classification)



library(dplyr)

summary_df <- cell_summary %>%
  mutate(
    sample_group = ident,
    cell_type = celltypes,
    group_label = paste0(sample_group, ifelse(cell_type == "Epithelial", "_Epi", "_All"))
  ) %>%
  group_by(group_label, classification) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group_label) %>%
  mutate(percent = count / sum(count) * 100)


library(ggplot2)

ggplot(summary_df, aes(x = group_label, y = count, fill = classification)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample Group", y = "Cell Count", fill = "Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Classification Counts by Sample")


epi_summary <- cell_summary %>%
  filter(celltypes == "Ductal") %>%
  group_by(ident, classification) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(ident) %>%
  mutate(percent = count / sum(count) * 100)

ggplot(epi_summary, aes(x = ident, y = percent, fill = classification)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Percentage of Epithelial Cells", fill = "Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Classification Profile of Epithelial Cells by Sample")



####################################################################



cna_signal <- vapply(
  genes_by_chr,
  function(g) {
    if (length(g) == 1) {
      outs[g, ]                       # keep as numeric vector
    } else {
      colMeans(outs[g, , drop = FALSE])
    }
  },
  numeric(ncol(outs))
)

cna_signal <- t(cna_signal)
chromosomes <- rownames(cna_signal)

ref_cells <- ref
ref_cells <- unlist(ref_cells, use.names = FALSE)

fit_params <- lapply(chromosomes, function(chr) {
  vals <- cna_signal[chr, ref_cells]
  fit  <- fitdistr(vals, "normal")
  c(mean = fit$estimate["mean"], sd = fit$estimate["sd"])
})
names(fit_params) <- chromosomes



zmat   <- matrix(NA_real_, nrow = nrow(cna_signal), ncol = ncol(cna_signal),
                 dimnames = dimnames(cna_signal))
p_mat  <- zmat
p_adj  <- zmat

for (chr in chromosomes) {
  mu <- fit_params[[chr]]["mean.mean"]
  sd <- fit_params[[chr]]["sd.sd"]
  z  <- (cna_signal[chr, ] - mu) / sd
  p  <- 2 * pnorm(-abs(z))                     # two-sided
  
  zmat[chr, ]  <- z
  p_mat[chr, ] <- p
  p_adj[chr, ] <- p.adjust(p, method = "BH")
}



call_matrix_chr <- matrix(NA_character_, nrow = nrow(cna_signal), 
                          ncol = ncol(cna_signal), dimnames = dimnames(cna_signal))

is_event <- p_adj < 0.05
is_susp  <- !is_event & (p_mat < 0.05)

call_matrix_chr[is_event & (cna_signal > 0)] <- "real_AMP"
call_matrix_chr[is_event & (cna_signal < 0)] <- "real_DEL"

call_matrix_chr[is_susp  & (cna_signal > 0)] <- "susp_AMP"
call_matrix_chr[is_susp  & (cna_signal < 0)] <- "susp_DEL"



cna_data <- read.delim("/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/matched/CNA_Genes.txt", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cna_data <- cna_data[cna_data$Profiled.Samples > 30, ]

panel_amp_genes <- intersect(cna_data$Gene[cna_data$CNA == "AMP"], rownames(outs))
panel_del_genes <- intersect(cna_data$Gene[cna_data$CNA == "HOMDEL"], rownames(outs))



cna_signal <- outs  # rows = genes, cols = cells
rownames(cna_signal) <- unlist(genes_by_chr, use.names = FALSE)

# Sort gene_order and cna_signal together
gene_order <- gene_order[match(rownames(cna_signal), gene_order$gene_id), ]
stopifnot(all.equal(rownames(cna_signal), gene_order$gene_id))

ref_cells <- unlist(ref, use.names = FALSE)
window_size <- 50

panel_genes <- unique(intersect(c(plane <- panel_amp_genes, panel_del_genes), rownames(cna_signal)))
amp_rows    <- intersect(panel_amp_genes, rownames(cna_signal))
del_rows    <- intersect(panel_del_genes, rownames(cna_signal))

# Small matrices: rows = panel genes, cols = all cells
zmat <- matrix(NA_real_, nrow = length(panel_genes), ncol = ncol(cna_signal),
               dimnames = list(panel_genes, colnames(cna_signal)))
p_mat <- zmat
p_adj <- zmat

# Precompute: map gene -> row index in cna_signal
gene_to_idx <- setNames(seq_len(nrow(cna_signal)), rownames(cna_signal))

# Sliding window by chromosome, but evaluate only panel genes
for (chr in unique(gene_order$chromosome)) {
  chr_gene_ids_all <- gene_order$gene_id[gene_order$chromosome == chr]
  chr_idx_all      <- gene_to_idx[intersect(chr_gene_ids_all, rownames(cna_signal))]
  chr_idx_all      <- as.integer(chr_idx_all[!is.na(chr_idx_all)])
  if (length(chr_idx_all) == 0) next
  
  # Panel genes on this chromosome
  chr_panel_genes <- intersect(chr_gene_ids_all, panel_genes)
  if (length(chr_panel_genes) == 0) next
  
  # For each target gene on this chr
  for (g in chr_panel_genes) {
    global_idx <- gene_to_idx[[g]]
    if (is.null(global_idx)) next
    
    rel_idx <- match(global_idx, chr_idx_all)
    if (is.na(rel_idx)) next
    
    half_w     <- floor(window_size / 2)
    start_idx  <- max(1, rel_idx - half_w)
    end_idx    <- min(length(chr_idx_all), rel_idx + half_w)
    window_idx <- chr_idx_all[start_idx:end_idx]
    
    ref_vals  <- cna_signal[window_idx, ref_cells, drop = FALSE]
    flat_vals <- as.numeric(ref_vals)
    
    fit <- tryCatch(
      fitdistr(flat_vals, "normal"),
      error = function(e) {
        m <- mean(flat_vals, na.rm = TRUE)
        s <- stats::sd(flat_vals, na.rm = TRUE)
        list(estimate = c(mean = m, sd = ifelse(is.finite(s) && s > 0, s, 1e-6)))
      }
    )
    
    mu <- unname(fit$estimate["mean"])
    sd <- unname(fit$estimate["sd"])
    if (!is.finite(sd) || sd <= 0) sd <- 1e-6
    
    gene_vals <- cna_signal[global_idx, ]
    
    z <- (gene_vals - mu) / sd
    p <- 2 * pnorm(-abs(z))  # two-sided
    
    zmat[g, ]  <- z
    p_mat[g, ] <- p
    p_adj[g, ] <- p.adjust(p, method = "BH")
  }
}

# --- Calls only for the relevant direction per panel set ---
call_matrix <- matrix(NA_character_, nrow = length(panel_genes), ncol = ncol(cna_signal),
                      dimnames = list(panel_genes, colnames(cna_signal)))

# AMP genes: positive direction
if (length(amp_rows)) {
  amp_Z   <- zmat[amp_rows, , drop = FALSE]
  amp_P   <- p_mat[amp_rows, , drop = FALSE]
  amp_Pad <- p_adj[amp_rows, , drop = FALSE]
  amp_X   <- cna_signal[amp_rows, , drop = FALSE]
  
  is_event_amp <- (amp_Pad < 0.05) & (amp_Z > 0)
  is_susp_amp  <- (amp_Pad >= 0.05 | is.na(amp_Pad)) & (amp_P < 0.05) & (amp_Z > 0)
  
  call_matrix[amp_rows, ][is_event_amp] <- "real_AMP"
  call_matrix[amp_rows, ][is_susp_amp]  <- "susp_AMP"
}

# DEL genes: negative direction
if (length(del_rows)) {
  del_Z   <- zmat[del_rows, , drop = FALSE]
  del_P   <- p_mat[del_rows, , drop = FALSE]
  del_Pad <- p_adj[del_rows, , drop = FALSE]
  del_X   <- cna_signal[del_rows, , drop = FALSE]
  
  is_event_del <- (del_Pad < 0.05) & (del_Z < 0)
  is_susp_del  <- (del_Pad >= 0.05 | is.na(del_Pad)) & (del_P < 0.05) & (del_Z < 0)
  
  call_matrix[del_rows, ][is_event_del] <- "real_DEL"
  call_matrix[del_rows, ][is_susp_del]  <- "susp_DEL"
}




celltype_init <- meta[[ct_col]][match(colnames(outs), rownames(meta))]
names(celltype_init) <- colnames(outs)

cell_status <- rep(NA_character_, length(celltype_init))
names(cell_status) <- names(celltype_init)

has_real <- apply(call_matrix, 2, function(col) any(grepl("^real_", col)))
has_suspicious <- apply(call_matrix, 2, function(col) any(grepl("^susp_", col)))

has_real_chr <- apply(call_matrix_chr, 2, function(col) any(grepl("^real_", col)))
has_suspicious_chr <- apply(call_matrix_chr, 2, function(col) any(grepl("^susp_", col)))


for (id in names(cell_status)) {
  
  if ( has_real[id] && has_real_chr[id] ) {
    cell_status[id] <- "malignant"
    
  } else if ( !has_real[id] && !has_suspicious[id] &&
              !has_real_chr[id] && !has_suspicious_chr[id] ) {
    cell_status[id] <- "nonmalignant"
    
  } else {
    cell_status[id] <- "unresolved"
  }
  
}



event_long <- as.data.frame(call_matrix) %>% 
  tibble::rownames_to_column("chromosome") %>% 
  pivot_longer(-chromosome, names_to = "cell_id", values_to = "call")

cell_summary <- data.frame(
  cell_id   = names(cell_status),
  ident = meta$orig.ident, 
  classification = cell_status,
  celltypes = meta[[ct_col]], 
  stringsAsFactors = FALSE
)

cell_status[cell_summary$classification == "malignant"]

table(cell_summary[cell_summary$ident == "H_post_T1_biopsy", ]$classification)
table(cell_summary[cell_summary$ident == "H_post_T1_biopsy" & cell_summary$celltypes == "epithelial", ]$classification)
table(cell_summary[cell_summary$ident != "H_post_T1_biopsy", ]$classification)



library(dplyr)

summary_df <- cell_summary %>%
  mutate(
    cell_type = celltypes,
    group_label = ifelse(cell_type == "epithelial", "patient_H_T1_biopsy", "Reference")
  ) %>%
  group_by(group_label, classification) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(group_label) %>%
  mutate(percent = count / sum(count) * 100)


library(ggplot2)

ggplot(summary_df, aes(x = group_label, y = count, fill = classification)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), 
            position = position_stack(vjust = 0.5), 
            fontface = "bold", 
            color = "black", 
            size = 3.5) +
  labs(x = "Sample Group", y = "Cell Count", fill = "Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15))
