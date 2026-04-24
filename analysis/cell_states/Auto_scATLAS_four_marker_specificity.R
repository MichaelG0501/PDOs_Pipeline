#!/usr/bin/env Rscript

####################################################
# Auto_scATLAS_four_marker_specificity.R
# Re-calculates metrics from Seurat objects for exact consistency
# Creates 4 pages:
# 1. scATLAS (4 target genes)
# 2. PDO (4 target genes)
# 3. scATLAS (Top 5 PDO Surface Markers for Basal and SMG)
# 4. PDO (Top 5 PDO Surface Markers for Basal and SMG)
####################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table)

# 1. Define base target genes
base_targets <- data.frame(
  gene = c("MUC13", "CEACAM5", "ROR1", "PTPRG"),
  target_state_group = c("Basal set", "Basal set", "SMG-like set", "SMG-like set"),
  stringsAsFactors = FALSE
)

# Canonical state orders
pdo_state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "3CA EMT and Protein maturation"
)

sc_state_order <- c(
  "Classic Proliferative",
  "Basal to Intestinal Metaplasia",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "Immune Infiltrating",
  "3CA EMT and Protein maturation"
)

# 2. Read Top 5 PDO SURFACE markers
# Path from AGENTS.md / discovery
surface_marker_path <- "PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_ranked.csv"
surface_markers <- fread(surface_marker_path)

# Extract top 5 for Basal and SMG (using PDO surface marker ranking)
pdo_surface_top5_basal <- surface_markers %>%
  filter(state == "Basal to Intest. Meta") %>%
  head(5) %>%
  pull(gene)

pdo_surface_top5_smg <- surface_markers %>%
  filter(state == "SMG-like Metaplasia") %>%
  head(5) %>%
  pull(gene)

# Define the targets for pages 3 & 4
surface_targets <- data.frame(
  gene = c(pdo_surface_top5_basal, pdo_surface_top5_smg),
  target_state_group = c(rep("Basal set", 5), rep("SMG-like set", 5)),
  stringsAsFactors = FALSE
)

# Union of all genes to compute
all_target_genes <- unique(c(base_targets$gene, surface_targets$gene))

####################
# PDO Calculation
####################
message("Loading PDO data...")
pdos_all <- readRDS("PDOs_outs/PDOs_merged.rds")
state_labels_pdo <- readRDS("PDOs_outs/Auto_PDO_final_states.rds")

# Subset to target genes and valid cells
valid_cells_pdo <- intersect(colnames(pdos_all), names(state_labels_pdo))
genes_in_pdo <- intersect(all_target_genes, rownames(pdos_all))
pdos_all <- pdos_all[genes_in_pdo, valid_cells_pdo]
pdos_all$state <- as.character(state_labels_pdo[valid_cells_pdo])

# Standardize state names for PDO
pdos_all$state <- gsub("Basal to Intest. Meta", "Basal to Intestinal Metaplasia", pdos_all$state)
pdos_all$state <- gsub("3CA_EMT_and_Protein_maturation", "3CA EMT and Protein maturation", pdos_all$state)

# Subset to the 5 finalized states
pdos_all <- subset(pdos_all, subset = state %in% pdo_state_order)
pdos_all <- NormalizeData(pdos_all, verbose = FALSE)

message("Computing PDO metrics...")
pdo_res <- list()
for (g in genes_in_pdo) {
  for (s in pdo_state_order) {
    cells_in_state <- colnames(pdos_all)[pdos_all$state == s]
    if (length(cells_in_state) == 0) next
    
    # Sample-aware metrics
    samples_in_state <- unique(pdos_all$orig.ident[cells_in_state])
    sample_pcts <- c()
    sample_means <- c()
    
    for (samp in samples_in_state) {
      samp_cells <- intersect(cells_in_state, colnames(pdos_all)[pdos_all$orig.ident == samp])
      if (length(samp_cells) < 10) next
      
      expr_vals <- GetAssayData(pdos_all, assay = "RNA", layer = "data")[g, samp_cells]
      sample_pcts <- c(sample_pcts, mean(expr_vals > 0))
      sample_means <- c(sample_means, mean(expr_vals))
    }
    
    if (length(sample_pcts) > 0) {
      pdo_res[[length(pdo_res) + 1]] <- data.frame(
        gene = g,
        state = s,
        dataset = "PDO",
        pct_expr = median(sample_pcts),
        avg_expr = median(sample_means),
        stringsAsFactors = FALSE
      )
    }
  }
}
pdo_data <- bind_rows(pdo_res)

rm(pdos_all)
invisible(gc())

####################
# scATLAS Calculation
####################
message("Loading scATLAS data...")
sc_path <- "/rds/general/user/sg3723/projects/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds"
sc_meta_path <- "/rds/general/user/sg3723/projects/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Auto_final_states.rds"

sc_obj <- readRDS(sc_path)
state_labels_sc <- readRDS(sc_meta_path)

# Ensure states are in the object
sc_obj$state <- as.character(state_labels_sc[match(colnames(sc_obj), names(state_labels_sc))])

# Standardize state names for scATLAS
sc_obj$state <- gsub("3CA_EMT_and_Protein_maturation", "3CA EMT and Protein maturation", sc_obj$state)

# Subset to target genes and valid states
genes_in_sc <- intersect(all_target_genes, rownames(sc_obj))
sc_obj <- sc_obj[genes_in_sc, ]
sc_obj <- subset(sc_obj, subset = state %in% sc_state_order)
sc_obj <- NormalizeData(sc_obj, verbose = FALSE)

message("Computing scATLAS metrics...")
sc_res <- list()
for (g in genes_in_sc) {
  for (s in sc_state_order) {
    cells_in_state <- colnames(sc_obj)[sc_obj$state == s]
    if (length(cells_in_state) == 0) next
    
    # Sample-aware metrics (using orig.ident)
    samples_in_state <- unique(sc_obj$orig.ident[cells_in_state])
    sample_pcts <- c()
    sample_means <- c()
    
    for (samp in samples_in_state) {
      samp_cells <- intersect(cells_in_state, colnames(sc_obj)[sc_obj$orig.ident == samp])
      if (length(samp_cells) < 10) next
      
      expr_vals <- GetAssayData(sc_obj, assay = "RNA", layer = "data")[g, samp_cells]
      sample_pcts <- c(sample_pcts, mean(expr_vals > 0))
      sample_means <- c(sample_means, mean(expr_vals))
    }
    
    if (length(sample_pcts) > 0) {
      sc_res[[length(sc_res) + 1]] <- data.frame(
        gene = g,
        state = s,
        dataset = "scATLAS",
        pct_expr = median(sample_pcts),
        avg_expr = median(sample_means),
        stringsAsFactors = FALSE
      )
    }
  }
}
sc_data <- bind_rows(sc_res)

rm(sc_obj)
invisible(gc())

####################
# Plotting function
####################
create_plot <- function(data, targets_df, dataset_name, state_order, title, subtitle) {
  # Merge to keep only the selected targets
  plot_data <- data %>%
    filter(dataset == dataset_name) %>%
    inner_join(targets_df, by = "gene")
  
  # Normalize expression per gene (within this plot's data)
  plot_data <- plot_data %>%
    group_by(gene) %>%
    mutate(
      norm_expr = if(sd(avg_expr, na.rm=T) > 0) {
        (avg_expr - mean(avg_expr, na.rm=T)) / sd(avg_expr, na.rm=T)
      } else {
        0
      }
    ) %>%
    ungroup()
  
  # Ensure factors for order
  plot_data$state <- factor(plot_data$state, levels = state_order)
  # Maintain order of genes as they appear in targets_df
  plot_data$gene <- factor(plot_data$gene, levels = rev(unique(targets_df$gene)))
  plot_data$target_state_group <- factor(plot_data$target_state_group, levels = c("Basal set", "SMG-like set"))

  p <- ggplot(plot_data, aes(x = state, y = gene)) +
    geom_point(aes(size = pct_expr, color = norm_expr)) +
    facet_grid(target_state_group ~ ., scales = "free_y", space = "free_y") +
    scale_size_continuous(range = c(2, 12), labels = scales::percent, limits = c(0, 1)) +
    scale_color_gradientn(
      colors = c("#1D4E89", "#F8F4EC", "#B22222"), 
      name = "Median log-Expr\n(Z-score)",
      limits = c(-1.5, 1.5),
      oob = scales::squish
    ) +
    theme_bw() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "",
      y = "",
      size = "% Cells Expressing"
    ) +
    theme(
      strip.background = element_rect(fill = "white", color = "white"),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", size = 10),
      axis.text.y = element_text(face = "bold", size = 11),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey90", fill = NA),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 16)
    )
  
  return(p)
}

####################
# Create PDF Pages
####################
out_file <- "PDOs_outs/Auto_scATLAS_four_marker_specificity.pdf"
pdf(out_file, width = 11, height = 7, useDingbats = FALSE)

# Page 1: scATLAS (4 base genes)
print(create_plot(
  data = sc_data,
  targets_df = base_targets,
  dataset_name = "scATLAS",
  state_order = sc_state_order,
  title = "Marker state-specificity: scATLAS Ref",
  subtitle = "4 Selected Markers"
))

# Page 2: PDO (4 base genes)
print(create_plot(
  data = pdo_data,
  targets_df = base_targets,
  dataset_name = "PDO",
  state_order = pdo_state_order,
  title = "Marker state-specificity: PDO",
  subtitle = "4 Selected Markers"
))

# Page 3: scATLAS (Top 5 PDO Surface Markers)
print(create_plot(
  data = sc_data,
  targets_df = surface_targets,
  dataset_name = "scATLAS",
  state_order = sc_state_order,
  title = "Marker state-specificity: scATLAS Ref",
  subtitle = "Top 5 PDO Surface Markers for Basal and SMG States"
))

# Page 4: PDO (Top 5 PDO Surface Markers)
print(create_plot(
  data = pdo_data,
  targets_df = surface_targets,
  dataset_name = "PDO",
  state_order = pdo_state_order,
  title = "Marker state-specificity: PDO",
  subtitle = "Top 5 PDO Surface Markers for Basal and SMG States"
))

dev.off()
message("Success. Output saved to: ", out_file)
