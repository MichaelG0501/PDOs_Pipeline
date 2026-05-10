####################
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(writexl)

# Set directories
setwd("/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline")
out_dir <- "PDOs_outs/Auto_untreated_comparison"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load data
message("Loading metadata and states...")
meta <- readRDS("PDOs_outs/PDOs_all_meta.rds")
states <- readRDS("PDOs_outs/Auto_PDO_final_states.rds")

# Ensure cell names match
message("Aligning cells...")
common_cells <- intersect(rownames(meta), names(states))
meta <- meta[common_cells, ]
meta$Final_State <- states[common_cells]

# Filter for relevant samples
samples_to_compare <- c("SUR1072_Untreated_PDO", "SUR1090_Untreated_PDO")
meta_sub <- meta %>% filter(orig.ident %in% samples_to_compare)

if (nrow(meta_sub) == 0) {
  stop("No cells found for the specified samples.")
}

# Define Canonical State Order
canonical_states <- c(
  "Classic Proliferative",
  "Basal to Intest. Meta",
  "SMG-like Metaplasia",
  "Stress-adaptive",
  "3CA_EMT_and_Protein_maturation",
  "Unresolved",
  "Hybrid"
)

# Colors from AGENTS.md + user preference (Hybrid = black)
state_colors <- c(
  "Classic Proliferative" = "#E41A1C",
  "Basal to Intest. Meta" = "#4DAF4A",
  "SMG-like Metaplasia" = "#FF7F00",
  "Stress-adaptive" = "#984EA3",
  "3CA_EMT_and_Protein_maturation" = "#377EB8",
  "Unresolved" = "#D9D9D9",
  "Hybrid" = "black"
)

# Calculate counts and proportions
message("Calculating proportions...")
summary_long <- meta_sub %>%
  group_by(orig.ident, Final_State) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(Proportion = Cell_Count / sum(Cell_Count)) %>%
  ungroup()

# Reshape for Excel
message("Preparing Excel output...")
# Ensure all states are present even if 0
all_states_df <- expand.grid(
  orig.ident = samples_to_compare,
  Final_State = canonical_states,
  stringsAsFactors = FALSE
)

summary_full <- all_states_df %>%
  left_join(summary_long, by = c("orig.ident", "Final_State")) %>%
  mutate(Cell_Count = replace_na(Cell_Count, 0),
         Proportion = replace_na(Proportion, 0))

summary_wide <- summary_full %>%
  mutate(State = factor(Final_State, levels = canonical_states)) %>%
  pivot_wider(
    id_cols = State,
    names_from = orig.ident,
    values_from = c(Proportion, Cell_Count)
  ) %>%
  arrange(State)

# Reorder columns as requested: SUR1072 prop, 1090 prop, 1072 count, 1090 count
summary_excel <- summary_wide %>%
  select(
    State,
    `SUR1072_Proportion` = Proportion_SUR1072_Untreated_PDO,
    `SUR1090_Proportion` = Proportion_SUR1090_Untreated_PDO,
    `SUR1072_Count` = Cell_Count_SUR1072_Untreated_PDO,
    `SUR1090_Count` = Cell_Count_SUR1090_Untreated_PDO
  )

# Save Excel
message("Saving Excel...")
write_xlsx(summary_excel, file.path(out_dir, "untreated_state_comparison.xlsx"))

# Plot Pie Charts
message("Generating pie charts...")
plot_data <- summary_long %>%
  mutate(Final_State = factor(Final_State, levels = rev(canonical_states)))

p_list <- lapply(samples_to_compare, function(s) {
  df <- plot_data %>% filter(orig.ident == s) %>% arrange(desc(Final_State))
  
  # Label placement: outside the pie
  df <- df %>%
    mutate(pos = cumsum(Proportion) - Proportion/2,
           Label = ifelse(Proportion > 0.005, paste0(round(Proportion * 100, 1), "%"), ""))
  
  ggplot(df, aes(x = 1, y = Proportion, fill = Final_State)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0.2) +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_manual(values = state_colors) +
    # Push labels out by setting x > 1.5
    geom_text(aes(x = 1.7, y = pos, label = Label), size = 3.5, fontface = "bold") +
    labs(title = gsub("_Untreated_PDO", " Untreated", s), fill = "Cell State") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right") +
    expand_limits(x = 1.8)
})

combined_plot <- p_list[[1]] + p_list[[2]] + plot_layout(guides = "collect")

message("Saving plots...")
ggsave(file.path(out_dir, "untreated_state_pie_charts.pdf"), combined_plot, width = 14, height = 7)
ggsave(file.path(out_dir, "untreated_state_pie_charts.png"), combined_plot, width = 14, height = 7)

message("Done. Outputs in: ", out_dir)
####################
