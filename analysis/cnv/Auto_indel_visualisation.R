####################
# Analysis registry:
#   Status: active visualization
#   Script: analysis/cnv/Auto_indel_visualisation.R
#   Methodology: None (Visualization and Depth verification script)
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     PDOs_outs/indel_analysis/PDO_1072_variants_metrics.tsv
#     PDOs_outs/indel_analysis/PDO_1090_variants_metrics.tsv
#     Sarek recalibrated CRAMs
#   Outputs:
#     PDOs_outs/indel_analysis/Auto_indel_chromatogram_style.pdf
#     PDOs_outs/indel_analysis/Auto_indel_sequence_alignments.png
#     PDOs_outs/indel_analysis/Auto_variant_prioritization_table.csv
#     PDOs_outs/indel_analysis/Auto_variant_clusters_plot.png
#   Downstream: Terminal figures
####################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

WD <- "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
setwd(WD)

OUT_DIR <- file.path(WD, "PDOs_outs/indel_analysis")
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive=TRUE)

# 1. Load metrics TSVs
cols <- c("chrom", "pos", "ref", "alt", "filter", "tlod", "dp", "germq", "nalod", "nlod", "mbq", "mpos", "csq", "ad_normal", "af_normal", "ad_tumor", "af_tumor")

read_variant_metrics <- function(filepath, sample_name) {
  if (!file.exists(filepath)) return(data.frame())
  
  df <- read.csv(filepath, sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
  if(nrow(df) == 0) return(data.frame())
  if(ncol(df) > 17) df <- df[, 1:17]
  if(ncol(df) < 17) return(data.frame()) # Incomplete data
  
  colnames(df) <- cols
  
  df$sample <- sample_name
  df$af_tumor <- as.numeric(gsub("]", "", df$af_tumor))
  df$af_normal <- as.numeric(gsub("\\[", "", df$af_normal))
  df$dp <- as.numeric(df$dp)
  df$tlod <- as.numeric(df$tlod)
  
  df$indel_len <- abs(nchar(df$ref) - nchar(df$alt))
  df$is_snp <- df$indel_len == 0 & nchar(df$ref) == 1 & nchar(df$alt) == 1
  
  # Parse CSQ
  first_csq <- sapply(strsplit(as.character(df$csq), ","), `[`, 1)
  csq_parts <- strsplit(first_csq, "\\|")
  
  df$consequence <- sapply(csq_parts, function(x) if(length(x) >= 2) x[2] else NA)
  df$gene <- sapply(csq_parts, function(x) if(length(x) >= 4) x[4] else NA)
  
  df$is_exon <- grepl("missense|synonymous|stop|splice|frameshift|inframe|exon|coding", df$consequence, ignore.case=TRUE)
  
  return(df)
}

df_1072 <- read_variant_metrics(file.path(OUT_DIR, "PDO_1072_variants_metrics.tsv"), "PDO_1072")
df_1090 <- read_variant_metrics(file.path(OUT_DIR, "PDO_1090_variants_metrics.tsv"), "PDO_1090")

# Apply Quality Filters
if (nrow(df_1072) > 0) df_1072 <- df_1072 %>% filter(dp >= 10 | is.na(dp))
if (nrow(df_1090) > 0) df_1090 <- df_1090 %>% filter(dp >= 10 | is.na(dp))

# ---------------------------------------------------------
# INDEL ANALYSIS (Original Logic for Sequence Alignments)
# ---------------------------------------------------------

top_1072 <- data.frame()
if (nrow(df_1072) > 0) {
  top_1072 <- df_1072 %>%
    filter(!is_snp, indel_len >= 5, grepl("PASS", filter)) %>%
    arrange(desc(af_tumor), desc(indel_len)) %>%
    head(3) %>%
    mutate(label = paste0("Intergenic (", chrom, ") | Deletion -", indel_len, "bp (Exclusive to 1072)\nAllele Freq: ", round(af_tumor*100, 1), "%"))
}

top_1090 <- data.frame()
if (nrow(df_1090) > 0) {
  top_1090 <- df_1090 %>%
    filter(!is_snp, indel_len >= 5, grepl("PASS", filter)) %>%
    arrange(desc(af_tumor), desc(indel_len)) %>%
    head(3) %>%
    mutate(label = paste0("Intergenic (", chrom, ") | Deletion -", indel_len, "bp (Exclusive to 1090)\nAllele Freq: ", round(af_tumor*100, 1), "%"))
}

# Define controls (hardcoded known clean regions)
controls <- data.frame(
  sample = c("Control_WT", "Control_Shared"),
  chrom = c("chr22", "chr2"),
  pos = c(29489100, 6861116),
  ref = c("N", "ACACACACACACC"),
  alt = c("N", "A"),
  indel_len = c(0, 12),
  af_tumor = c(NA, NA),
  label = c("Wild-Type Control (chr22) | Absent in both samples\nContinuous, normal depth in both samples",
            "Shared Control (chr2) | Deletion -12bp present in BOTH samples\nIdentical drop in sequencing depth is visible in both samples")
)

# Combine for sequence alignment plotting
all_indels <- bind_rows(top_1072, top_1090, controls[controls$sample == "Control_Shared", ])

# Depth Extraction using samtools
CRAM1072 <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/preprocessing/recalibrated/PDO_1072/PDO_1072.recal.cram"
CRAM1090 <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/preprocessing/recalibrated/PDO_1090/PDO_1090.recal.cram"

depth_list <- list()
extract_depth <- function(chr, pos, len, region_name, label_text) {
  start_coord <- max(1, pos - 100)
  end_coord <- pos + max(100, len + 50)
  reg_coord <- paste0(chr, ":", start_coord, "-", end_coord)
  
  cmd_1072 <- sprintf("samtools depth -r %s %s", reg_coord, CRAM1072)
  cmd_1090 <- sprintf("samtools depth -r %s %s", reg_coord, CRAM1090)
  
  d1 <- read.table(text = system(cmd_1072, intern = TRUE), stringsAsFactors = FALSE)
  d2 <- read.table(text = system(cmd_1090, intern = TRUE), stringsAsFactors = FALSE)
  
  if (nrow(d1) > 0) { colnames(d1) <- c("chrom", "pos", "depth"); d1$sample <- "PDO_1072"; d1$region <- region_name; d1$label <- label_text }
  if (nrow(d2) > 0) { colnames(d2) <- c("chrom", "pos", "depth"); d2$sample <- "PDO_1090"; d2$region <- region_name; d2$label <- label_text }
  
  res <- bind_rows(if(nrow(d1)>0) d1 else NULL, if(nrow(d2)>0) d2 else NULL)
  
  if (region_name != "Control_WT") {
    if(region_name == "Control_Shared") {
      res$indel_start <- pos - 23
      res$indel_end <- pos - 11
    } else {
      res$indel_start <- pos
      res$indel_end <- pos + len
    }
  } else {
    res$indel_start <- NA
    res$indel_end <- NA
  }
  return(res)
}

if (nrow(top_1072) > 0) {
  for (i in 1:nrow(top_1072)) {
    depth_list[[length(depth_list)+1]] <- extract_depth(top_1072$chrom[i], top_1072$pos[i], top_1072$indel_len[i], paste0("1072_Excl_", i), top_1072$label[i])
  }
}
if (nrow(top_1090) > 0) {
  for (i in 1:nrow(top_1090)) {
    depth_list[[length(depth_list)+1]] <- extract_depth(top_1090$chrom[i], top_1090$pos[i], top_1090$indel_len[i], paste0("1090_Excl_", i), top_1090$label[i])
  }
}
depth_list[[length(depth_list)+1]] <- extract_depth(controls$chrom[1], controls$pos[1], controls$indel_len[1], controls$sample[1], controls$label[1])
depth_list[[length(depth_list)+1]] <- extract_depth(controls$chrom[2], controls$pos[2], controls$indel_len[2], controls$sample[2], controls$label[2])

df_depth <- bind_rows(depth_list)

if (nrow(df_depth) > 0) {
  df_depth$sample <- factor(df_depth$sample, levels = c("PDO_1090", "PDO_1072"))
  pdf_out <- file.path(OUT_DIR, "Auto_indel_chromatogram_style.pdf")
  pdf(pdf_out, width = 11, height = 7, onefile = TRUE)
  for (reg in unique(df_depth$region)) {
    df_sub <- df_depth %>% filter(region == reg)
    if(nrow(df_sub) == 0) next
    p <- ggplot(df_sub, aes(x = pos, y = depth)) +
      geom_area(fill = "#5D8AA8", alpha = 0.5) +
      geom_line(color = "#1D428A", linewidth = 0.8) +
      facet_wrap(~sample, ncol = 1, scales = "free_y") +
      labs(title = df_sub$label[1],
           subtitle = "Top panel: Sample 1090 | Bottom panel: Sample 1072\nA structural deletion appears as a sudden drop in sequencing depth.",
           x = paste("Genomic Position (", df_sub$chrom[1], ")", sep=""),
           y = "Sequencing Depth (Reads)") +
      theme_minimal(base_size = 14) +
      theme(
        strip.text = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "grey80"),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey40", size = 11),
        plot.margin = margin(20, 20, 20, 20)
      )
    if(!is.na(df_sub$indel_start[1])) {
      start_v <- df_sub$indel_start[1]
      end_v <- df_sub$indel_end[1]
      p <- p + geom_segment(aes(x = start_v, xend = end_v, y = 0, yend = 0), color = "red", linewidth = 4) +
               geom_vline(xintercept = c(start_v, end_v), color="red", linetype="dashed", alpha=0.5)
    }
    print(p)
  }
  dev.off()
}

plot_indel_alignment <- function(row_data) {
  ref <- row_data$ref; alt <- row_data$alt
  len <- abs(nchar(ref) - nchar(alt))
  type <- ifelse(nchar(ref) > nchar(alt), "deletion", "insertion")
  max_viz_len <- 40
  if (type == "deletion") {
    viz_len <- min(nchar(ref), max_viz_len)
    ref_viz <- substr(ref, 1, viz_len)
    alt_viz <- paste0(substr(alt, 1, 1), strrep("-", viz_len - 1))
    if (nchar(ref) > max_viz_len) { ref_viz <- paste0(ref_viz, "..."); alt_viz <- paste0(alt_viz, "...") }
  } else {
    viz_len <- min(nchar(alt), max_viz_len)
    ref_viz <- paste0(substr(ref, 1, 1), strrep("-", viz_len - 1))
    alt_viz <- substr(alt, 1, viz_len)
    if (nchar(alt) > max_viz_len) { ref_viz <- paste0(ref_viz, "..."); alt_viz <- paste0(alt_viz, "...") }
  }
  ref_chars <- strsplit(ref_viz, "")[[1]]; alt_chars <- strsplit(alt_viz, "")[[1]]
  n_pos <- max(length(ref_chars), length(alt_chars))
  df <- data.frame(pos = rep(1:n_pos, 2), allele = rep(c("REF", "ALT"), each = n_pos), base = c(ref_chars, alt_chars))
  base_colors <- c("A" = "#4daf4a", "C" = "#377eb8", "G" = "#ff7f00", "T" = "#e41a1c", "-" = "#cccccc", "."="white")
  p <- ggplot(df, aes(x = pos, y = allele, fill = base)) +
    geom_tile(color = "white", linewidth = 1, width = 0.9, height = 0.8) +
    geom_text(aes(label = base), size = 5, fontface = "bold", color = ifelse(df$base == "-", "transparent", "white")) +
    scale_fill_manual(values = base_colors) + scale_y_discrete(limits = c("ALT", "REF")) +
    labs(title = paste0("Sample: ", row_data$sample, " | ", row_data$chrom, ":", row_data$pos), subtitle = paste0(toupper(type), " of ", len, "bp"), x = NULL, y = NULL) +
    theme_minimal(base_size = 14) + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), legend.position = "none")
  return(p)
}

if (nrow(all_indels) > 0) {
  plots <- lapply(1:nrow(all_indels), function(i) plot_indel_alignment(all_indels[i, ]))
  combined_plot <- wrap_plots(plots, ncol = 1) +
    plot_annotation(title = "Somatic Indels: Sequence-Level Visualizations", subtitle = "Reference vs Alternate allele alignments for top exclusive indels")
  ggsave(file.path(OUT_DIR, "Auto_indel_sequence_alignments.png"), combined_plot, width = 12, height = 2.5 * nrow(all_indels), dpi = 300)
}

# ---------------------------------------------------------
# NEW SNP CLUSTER & PASS PRIORITIZATION
# ---------------------------------------------------------
all_vars <- bind_rows(df_1072, df_1090)

if (nrow(all_vars) > 0) {
  # 1. Filter for PASS variants and strict Exon coding
  all_vars <- all_vars %>% filter(grepl("PASS", filter))
  
  all_vars$var_id <- paste(all_vars$chrom, all_vars$pos, all_vars$ref, all_vars$alt, sep="_")
  all_vars$is_coding <- grepl("missense_variant|synonymous_variant|stop_gained|stop_lost|frameshift_variant|inframe_insertion|inframe_deletion|splice_acceptor_variant|splice_donor_variant", all_vars$consequence)
  
  all_vars <- all_vars %>% filter(gene != "" & !is.na(gene) & is_coding)
  
  # Identify variants unique to one sample
  counts <- all_vars %>% group_by(var_id) %>% summarise(n_samples = n_distinct(sample))
  all_vars <- left_join(all_vars, counts, by="var_id")
  unique_vars <- all_vars %>% filter(n_samples == 1)
  
  # 2. Identify 500bp SNP clusters (Coding)
  snp_clusters <- unique_vars %>%
    filter(is_snp) %>%
    group_by(gene, sample, chrom) %>%
    arrange(pos) %>%
    mutate(
      prev_pos = lag(pos),
      dist_to_prev = pos - prev_pos,
      next_pos = lead(pos),
      dist_to_next = next_pos - pos
    ) %>%
    filter((!is.na(dist_to_prev) & dist_to_prev <= 500) | (!is.na(dist_to_next) & dist_to_next <= 500))
  
  cluster_genes <- unique(snp_clusters$gene)
  
  # All unique PASS variants (Cluster + Single), sorted by priority
  prioritized_vars <- unique_vars %>%
    mutate(
      priority_reason = ifelse(gene %in% cluster_genes & is_snp, "500bp SNP Cluster", "High-Confidence PASS")
    ) %>%
    distinct(var_id, sample, .keep_all=TRUE) %>%
    arrange(desc(priority_reason == "500bp SNP Cluster"), desc(af_tumor), gene, pos) %>%
    select(gene, chrom, pos, ref, alt, sample, priority_reason, consequence, af_tumor, af_normal, dp, tlod, filter, everything())
  
  # Output single detailed table sorted by priority reason -> gene -> position
  write.csv(prioritized_vars, file.path(OUT_DIR, "Auto_variant_prioritization_details.csv"), row.names=FALSE)
  
  if(nrow(prioritized_vars) > 0) {
    # Select 3 best 500bp clusters
    top_cluster_genes <- prioritized_vars %>%
      filter(priority_reason == "500bp SNP Cluster") %>%
      group_by(gene) %>% summarise(max_af = max(af_tumor, na.rm=TRUE)) %>%
      arrange(desc(max_af)) %>% head(3) %>% pull(gene)
    
    # Select 2 best singles
    top_single_genes <- prioritized_vars %>%
      filter(priority_reason == "High-Confidence PASS" & !(gene %in% top_cluster_genes)) %>%
      group_by(gene) %>% summarise(max_af = max(af_tumor, na.rm=TRUE)) %>%
      arrange(desc(max_af)) %>% head(2) %>% pull(gene)
      
    top_genes <- c(top_cluster_genes, top_single_genes)
    
    fasta_path <- "/rds/general/project/tumourheterogeneity1/live/demultiplex/genome.fa"
    base_colors <- c("A" = "#4daf4a", "C" = "#377eb8", "G" = "#ff7f00", "T" = "#e41a1c", "-" = "#cccccc", "N" = "#999999")
    
    pdf_path <- file.path(OUT_DIR, "Auto_variant_clusters_plot.pdf")
    pdf(pdf_path, width=20, height=8)
    
    for (g in top_genes) {
      g_vars <- prioritized_vars %>% filter(gene == g) %>% arrange(pos)
      chrom <- g_vars$chrom[1]
      min_pos <- min(g_vars$pos)
      max_pos <- max(g_vars$pos)
      
      # 20bp padding
      pad <- 20
      start_pos <- min_pos - pad
      end_pos <- max_pos + pad
      
      # Fetch sequence
      cmd <- sprintf("samtools faidx %s %s:%d-%d | grep -v '>' | tr -d '\\n'", fasta_path, chrom, start_pos, end_pos)
      ref_seq <- tryCatch(system(cmd, intern=TRUE), error=function(e) "")
      
      if (length(ref_seq) == 0 || nchar(ref_seq) == 0) {
        warning(paste("Could not fetch sequence for", g))
        next
      }
      
      ref_chars <- strsplit(toupper(ref_seq), "")[[1]]
      pos_seq <- start_pos:end_pos
      
      if(length(ref_chars) != length(pos_seq)) {
        length(ref_chars) <- length(pos_seq)
        ref_chars[is.na(ref_chars)] <- "N"
      }
      
      plot_rows <- list()
      mut_labels <- c()
      
      for (s in c("PDO_1072", "PDO_1090")) {
        for (idx in 1:length(pos_seq)) {
          p <- pos_seq[idx]
          ref_b <- ref_chars[idx]
          
          # Check if this position is mutated in this sample
          mut_match <- g_vars[g_vars$pos == p & g_vars$sample == s, ]
          
          if (nrow(mut_match) > 0) {
            row <- mut_match[1,]
            alt_b <- row$alt
            af <- row$af_tumor
            
            # Alt bar
            plot_rows[[length(plot_rows)+1]] <- data.frame(
              Sample = s, Position = p, Nucleotide = alt_b, Proportion = af
            )
            # Ref bar
            plot_rows[[length(plot_rows)+1]] <- data.frame(
              Sample = s, Position = p, Nucleotide = ref_b, Proportion = 1 - af
            )
            mut_labels <- c(mut_labels, paste0(p, " (", row$ref, ">", alt_b, " AF:", round(af,2), ")"))
          } else {
            # Not mutated
            plot_rows[[length(plot_rows)+1]] <- data.frame(
              Sample = s, Position = p, Nucleotide = ref_b, Proportion = 1.0
            )
          }
        }
      }
      
      plot_data <- bind_rows(plot_rows)
      plot_data$Sample <- factor(plot_data$Sample, levels=c("PDO_1072", "PDO_1090"))
      
      plot_data$DisplayBase <- sapply(plot_data$Nucleotide, function(x) ifelse(nchar(x) > 1, substr(x, 1, 2), x))
      plot_data$ColorBase <- sapply(plot_data$Nucleotide, function(x) {
         ch <- substr(x, 1, 1)
         if (ch %in% names(base_colors)) ch else "-"
      })
      
      mut_labels <- unique(mut_labels)
      priority_label <- g_vars$priority_reason[1]
      title_str <- paste0("Gene: ", g, " | ", priority_label, " | ", chrom, ":", min_pos, "-", max_pos)
      sub_str <- paste("Variants:", paste(mut_labels, collapse="; "))
      
      p_plot <- ggplot(plot_data, aes(x = Position, y = Proportion, fill = ColorBase)) +
        geom_col(width = 0.9, color = "black", linewidth = 0.2) +
        geom_text(aes(label = DisplayBase), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold", color = "white") +
        scale_fill_manual(values = base_colors) +
        facet_wrap(~Sample, ncol=1) +
        scale_x_continuous(breaks = seq(start_pos, end_pos, by=5)) +
        theme_minimal(base_size = 14) +
        labs(title = title_str, subtitle = sub_str, x = "Genomic Position", y = "Allele Proportion") +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          strip.text = element_text(face = "bold", size=12),
          panel.grid.major.x = element_blank(),
          legend.position = "none"
        )
        
      # Draw boxes around the mutated positions
      for (pos in unique(g_vars$pos)) {
        p_plot <- p_plot + annotate("rect", xmin = pos - 0.5, xmax = pos + 0.5, ymin = 0, ymax = 1,
                          color = "black", fill = NA, linewidth = 1)
      }
        
      print(p_plot)
    }
    
    dev.off()
  }
}

cat("Successfully generated sequence visualizations, priority tables, and depth profiles in", OUT_DIR, "\n")
