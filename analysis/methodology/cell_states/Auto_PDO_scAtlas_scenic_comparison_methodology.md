# Auto_PDO_scAtlas_scenic_comparison Methodology

## Overview
This document outlines the methodology for comparing SCENIC regulon activities between the single-cell Atlas (scRef) and Patient-Derived Organoids (PDOs). The primary goal is to assess whether the gene regulatory networks (GRNs) underpinning the canonical epithelial states are conserved ex vivo.

## Core Statistical Units

1. **Raw AUCell (AUC)**: The fundamental output of SCENIC. It evaluates the recovery of a regulon's target genes within the expression ranking of each cell. AUC is an absolute metric that is broadly comparable across datasets.
2. **Regulon Specificity Score (RSS)**: Calculated based on the Jensen-Shannon divergence of AUCell scores across cell states. RSS measures how specific a regulon's activity is to a given state.
3. **Specificity Gap**: The primary ranking metric used in this workflow. It is defined as:
   `Gap = RSS(State of Interest) - max(RSS(All Other States))`
   This strictly quantifies the exclusivity of a regulon to a specific state. A positive gap indicates the regulon is uniquely active in that state compared to all others.

## Selection and Plotting Thresholds

### 1. Broad Exploratory Heatmaps
**Files**: `Auto_scenic_comparison_MP_heatmap.pdf`, `Auto_scenic_comparison_State_heatmap.pdf`
- **Selection**: No strict filtering is applied. These heatmaps plot the intersection of **all common regulons** present in both datasets.
- **Scaling Method**: To visually align both datasets on the same relative scale and emphasize structural state-to-state similarities, the combined RSS matrices undergo **row-wise Z-scoring** (`t(scale(t(mat)))`). 
- **Important Visualization Limitation**: Z-scoring shifts the mean of every regulon to 0 and its standard deviation to 1. **This severely amplifies weak/inactive regulons**. For example, if a regulon like `TP73` has universally negligible absolute RSS across all states, but is fractionally higher in the Proliferative state, the Z-score formula will still force it to appear as a high-intensity peak (bright red). **Always cross-reference the Excel workbook for absolute magnitude to verify if a "hot" regulon on this heatmap is actually biologically active.**

### 2. Top-Tier Highlights
**Files**: `Auto_scenic_comparison_RSS_heatmaps.pdf` (Multi-page PDF), `Auto_scRef_PDO_scenic_top5_markers.xlsx`
- **Selection**: Explicitly filtered to the **top 5 regulons per state**. 
- **Ranking**: Ranked exclusively by the **Specificity Gap**. For the combined cross-dataset visualization, regulons must have a positive gap in *both* datasets (`sc_vals > 0 & pdo_vals > 0`), and are ranked by the mean of their gaps.
- **Scaling Method**: These heatmaps plot the **raw, unscaled Specificity Gap**, bypassing the Z-score amplification issue.

## State Mapping
To account for dataset-specific nomenclature, states are mapped logically based on their hierarchical tier:
- `Classic proliferation` (Shared)
- `Stress-adaptive` (Shared)
- `Cancer-cell immune mimicry` (scAtlas) maps independently (no strict ex vivo equivalent).
- `Squamous-to-intestinal` (scAtlas) aligns with `Columnar-to-intestinal` (PDOs) under a shared `Intestinal Metaplasia` envelope.
- `Glandular-to-intestinal` (scAtlas) aligns with `Glandular differentiation` (PDOs). 

## Storage and Caching
- **Canonical Script**: `analysis/cell_states/Auto_PDO_scAtlas_scenic_comparison.R`
- **Inputs**: Requires the complete SCENIC outputs from both the live `scRef_Pipeline` and `PDOs_Pipeline` `final_mp_scenic` runs.
- **Outputs**: All outputs must be saved to the persistent `PDOs_outs/final_mp_scenic` live directory. No outputs from this comparison are currently used downstream by other analytical scripts.
