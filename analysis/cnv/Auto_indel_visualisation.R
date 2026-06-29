####################
# Analysis registry:
#   Status: active visualization
#   Script: analysis/cnv/Auto_indel_visualisation.R
#   Methodology: None (Visualization and Depth verification script)
#   Map: analysis/ANALYSIS_MAP.md
#   Inputs:
#     PDOs_outs/indel_analysis/PDO_1072_indels_metrics.tsv
#     PDOs_outs/indel_analysis/PDO_1090_indels_metrics.tsv
#     Sarek recalibrated CRAMs
#   Outputs:
#     PDOs_outs/indel_analysis/Auto_indel_chromatogram_style.pdf
#     PDOs_outs/indel_analysis/Auto_indel_sequence_alignments.png
#   Downstream: Terminal figures
####################

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

WD <- "/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
setwd(WD)

OUT_DIR <- file.path(WD, "PDOs_outs/indel_analysis")
if(!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive=TRUE)

# 1. Load metrics TSVs
cols <- c("chrom", "pos", "ref", "alt", "filter", "tlod", "dp", "germq", "nalod", "nlod", "mbq", "mpos", "ad_normal", "af_normal", "ad_tumor", "af_tumor")
df_1072 <- read.csv(file.path(OUT_DIR, "PDO_1072_indels_metrics.tsv"), sep="\t", header=FALSE, fill=TRUE)
df_1090 <- read.csv(file.path(OUT_DIR, "PDO_1090_indels_metrics.tsv"), sep="\t", header=FALSE, fill=TRUE)
if(ncol(df_1072) > 16) df_1072 <- df_1072[, 1:16]
if(ncol(df_1090) > 16) df_1090 <- df_1090[, 1:16]
colnames(df_1072) <- cols
colnames(df_1090) <- cols

# Clean columns
df_1072$af_tumor <- as.numeric(gsub("]", "", df_1072$af_tumor))
df_1090$af_tumor <- as.numeric(gsub("]", "", df_1090$af_tumor))

df_1072$indel_len <- abs(nchar(df_1072$ref) - nchar(df_1072$alt))
df_1090$indel_len <- abs(nchar(df_1090$ref) - nchar(df_1090$alt))

# 2. Select top 3 large indels (>5bp) exclusive to each
top_1072 <- df_1072 %>%
  filter(indel_len >= 5) %>%
  arrange(desc(af_tumor), desc(indel_len)) %>%
  head(3) %>%
  mutate(sample = "PDO_1072", label = paste0("Intergenic (", chrom, ") | Deletion -", indel_len, "bp (Exclusive to 1072)\nAllele Freq: ", round(af_tumor*100, 1), "%"))

top_1090 <- df_1090 %>%
  filter(indel_len >= 5) %>%
  arrange(desc(af_tumor), desc(indel_len)) %>%
  head(3) %>%
  mutate(sample = "PDO_1090", label = paste0("Intergenic (", chrom, ") | Deletion -", indel_len, "bp (Exclusive to 1090)\nAllele Freq: ", round(af_tumor*100, 1), "%"))

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

# 3. Extract depth using samtools
CRAM1072 <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/preprocessing/recalibrated/PDO_1072/PDO_1072.recal.cram"
CRAM1090 <- "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/preprocessing/recalibrated/PDO_1090/PDO_1090.recal.cram"

depth_list <- list()

extract_depth <- function(chr, pos, len, region_name, label_text) {
  start_coord <- max(1, pos - 100)
  end_coord <- pos + max(100, len + 50)
  reg_coord <- paste0(chr, ":", start_coord, "-", end_coord)
  
  cat("Extracting depth for", region_name, "at", reg_coord, "...\n")
  
  cmd_1072 <- sprintf("samtools depth -r %s %s", reg_coord, CRAM1072)
  cmd_1090 <- sprintf("samtools depth -r %s %s", reg_coord, CRAM1090)
  
  d1 <- read.table(text = system(cmd_1072, intern = TRUE), stringsAsFactors = FALSE)
  d2 <- read.table(text = system(cmd_1090, intern = TRUE), stringsAsFactors = FALSE)
  
  if (nrow(d1) > 0) { colnames(d1) <- c("chrom", "pos", "depth"); d1$sample <- "PDO_1072"; d1$region <- region_name; d1$label <- label_text }
  if (nrow(d2) > 0) { colnames(d2) <- c("chrom", "pos", "depth"); d2$sample <- "PDO_1090"; d2$region <- region_name; d2$label <- label_text }
  
  res <- bind_rows(if(nrow(d1)>0) d1 else NULL, if(nrow(d2)>0) d2 else NULL)
  
  # Add highlight coordinates
  if (region_name != "Control_WT") {
    # Adjust for left-alignment by aligner in CA repeat
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

# Run extraction
for (i in 1:nrow(top_1072)) {
  depth_list[[length(depth_list)+1]] <- extract_depth(top_1072$chrom[i], top_1072$pos[i], top_1072$indel_len[i], paste0("1072_Excl_", i), top_1072$label[i])
}
for (i in 1:nrow(top_1090)) {
  depth_list[[length(depth_list)+1]] <- extract_depth(top_1090$chrom[i], top_1090$pos[i], top_1090$indel_len[i], paste0("1090_Excl_", i), top_1090$label[i])
}
depth_list[[length(depth_list)+1]] <- extract_depth(controls$chrom[1], controls$pos[1], controls$indel_len[1], controls$sample[1], controls$label[1])
depth_list[[length(depth_list)+1]] <- extract_depth(controls$chrom[2], controls$pos[2], controls$indel_len[2], controls$sample[2], controls$label[2])

df_depth <- bind_rows(depth_list)

# 4. Plot Depth Coverage
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
    p <- p + geom_segment(aes(x = start_v, xend = end_v, y = 0, yend = 0),
                          color = "red", linewidth = 4) +
             geom_vline(xintercept = c(start_v, end_v), color="red", linetype="dashed", alpha=0.5)
  }
  
  print(p)
}
dev.off()

# 5. Plot Sequence Alignments
plot_indel_alignment <- function(row_data) {
  ref <- row_data$ref
  alt <- row_data$alt
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
  
  ref_chars <- strsplit(ref_viz, "")[[1]]
  alt_chars <- strsplit(alt_viz, "")[[1]]
  n_pos <- max(length(ref_chars), length(alt_chars))
  
  df <- data.frame(
    pos = rep(1:n_pos, 2),
    allele = rep(c("REF", "ALT"), each = n_pos),
    base = c(ref_chars, alt_chars)
  )
  
  df$color_type <- "Match"
  if (type == "deletion") {
    df$color_type[df$allele == "REF" & df$pos > 1 & df$pos <= len + 1] <- "Deleted"
    df$color_type[df$allele == "ALT" & df$base == "-"] <- "Gap"
  } else {
    df$color_type[df$allele == "ALT" & df$pos > 1 & df$pos <= len + 1] <- "Inserted"
    df$color_type[df$allele == "REF" & df$base == "-"] <- "Gap"
  }
  
  base_colors <- c("A" = "#4daf4a", "C" = "#377eb8", "G" = "#ff7f00", "T" = "#e41a1c", "-" = "#cccccc", "."="white")
  
  title_text <- paste0("Sample: ", row_data$sample, " | ", row_data$chrom, ":", row_data$pos)
  subtitle_text <- paste0(toupper(type), " of ", len, "bp")
  
  p <- ggplot(df, aes(x = pos, y = allele, fill = base)) +
    geom_tile(color = "white", linewidth = 1, width = 0.9, height = 0.8) +
    geom_text(aes(label = base), size = 5, fontface = "bold", color = ifelse(df$base == "-", "transparent", "white")) +
    scale_fill_manual(values = base_colors) +
    scale_y_discrete(limits = c("ALT", "REF")) +
    labs(title = title_text, subtitle = subtitle_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )
    
  return(p)
}

plots <- lapply(1:nrow(all_indels), function(i) plot_indel_alignment(all_indels[i, ]))
combined_plot <- wrap_plots(plots, ncol = 1) +
  plot_annotation(
    title = "Somatic Indels: Sequence-Level Visualizations",
    subtitle = "Reference vs Alternate allele alignments for top exclusive indels",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey40")
    )
  )

png_out <- file.path(OUT_DIR, "Auto_indel_sequence_alignments.png")
ggsave(png_out, combined_plot, width = 12, height = 2.5 * nrow(all_indels), dpi = 300)

cat("Successfully generated sequence visualizations and depth profiles in", OUT_DIR, "\n")
