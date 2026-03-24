####################
# Auto_pdo_overall_state_proportions.R
# Overall proportion barplot for PDO states (Approach B).
#
# Input:
#   PDOs_outs/PDOs_merged.rds
#   PDOs_outs/Auto_PDO_final_states.rds
#
# Output:
#   PDOs_outs/Auto_pdo_overall_state_proportions.pdf
####################

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# setwd to PDOs_outs
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs")

# Load data
message("Loading data ...")
tmdata_all <- readRDS("PDOs_merged.rds")
state_B <- readRDS("Auto_PDO_final_states.rds")

# Common cells
common_cells <- intersect(names(state_B), Cells(tmdata_all))
state_B <- state_B[common_cells]

# State order and colors (scRef nomenclature - PDO data uses abbreviated names)
group_cols <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "Stress-adaptive"       = "#984EA3",
  "SMG-like Metaplasia"   = "#FF7F00",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved"            = "grey80",
  "Hybrid"                = "black"
)

# Ensure correct order for visualization
state_order <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "Stress-adaptive",
  "SMG-like Metaplasia",
  "3CA_EMT_and_Protein_maturation",
  "Unresolved",
  "Hybrid"
)

# Calculate proportions
prop_df <- data.frame(
  state = factor(as.character(state_B), levels = state_order),
  stringsAsFactors = FALSE
)

overall <- prop_df %>%
  count(state, .drop = FALSE) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  mutate(label = sprintf("%.1f%%", pct))

# Plotting overall barplot (single bar)
message("Generating overall barplot ...")
p <- ggplot(overall, aes(x = "Overall", y = pct, fill = state)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), 
            size = 6, 
            fontface = "bold",
            color = ifelse(overall$state %in% c("Hybrid"), "white", "black")) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  labs(title = "PDO State Proportions",
       subtitle = paste0("Total cells: ", sum(overall$n)),
       x = NULL, y = "% of cells", fill = "State") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5), 
        legend.title = element_text(size = 16),
        legend.text  = element_text(size = 14)) +
  coord_cartesian(expand = FALSE)

# Save output
pdf("Auto_pdo_overall_state_proportions.pdf", width = 7, height = 9)
print(p)
dev.off()

message("Finished. Plot saved to PDOs_outs/Auto_pdo_overall_state_proportions.pdf")
