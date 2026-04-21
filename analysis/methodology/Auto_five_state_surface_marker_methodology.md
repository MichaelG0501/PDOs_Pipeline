# Five-State PDO Surface Marker Methodology

This document describes how `analysis/cell_states/Auto_five_state_surface_markers.R` prioritizes PDO marker genes for FACS-oriented cell-surface targeting.

It is written against the current code and the current upstream marker outputs from:

- `analysis/cell_states/Auto_five_state_markers.R`
- `analysis/cell_states/Auto_five_state_surface_markers.R`

The description below is operational. It documents the implemented logic, not a simplified conceptual version.

## 1. Marker Inputs Reused From The PDO Marker Workflow

### 1.1 Reused marker tables

The script does not recompute PDO differential expression.

It reuses:

- `PDOs_outs/Auto_five_state_markers/Auto_five_state_marker_summary.csv`
- `PDOs_outs/Auto_five_state_markers/Auto_five_state_markers_ranked.csv`

Each row is a gene-state pair from the five finalized PDO states.

### 1.2 Initial candidate reuse filter

Before any surface annotation is applied, the script keeps only gene-state pairs that already passed the upstream PDO marker logic in a minimal way:

- `hit_sample_n > 0`
- `best_state_match == TRUE`
- `specificity_gap > 0`

This means the gene:

- was detected as an upregulated hit in at least one eligible sample for that state
- had its highest median sample-level expression in that state
- had a positive target-vs-best-off-state sample-level expression gap

This keeps the surface-marker workflow anchored to the existing PDO-only marker definition instead of starting from all genes.

## 2. PDO Expression Object Used For FACS Metrics

### 2.1 Input objects

The script loads:

- `PDOs_outs/PDOs_merged.rds`
- `PDOs_outs/Auto_PDO_final_states.rds`

It restricts the analysis to the five finalized PDO states:

- `Classic Proliferative`
- `Basal to Intest. Meta`
- `Stress-adaptive`
- `SMG-like Metaplasia`
- `3CA_EMT_and_Protein_maturation`

### 2.2 Expression layers used

The script reads the `RNA` assay:

- `counts` layer for fraction of cells expressing a gene
- `data` layer for mean normalized log-expression

No extra normalization is performed in this surface-marker step. It uses the expression representation already stored in `PDOs_merged.rds`.

## 3. External Annotation References

### 3.1 Download and cache behavior

If the external annotation files are missing, the script downloads and caches them.

Default cache directory:

- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_pdo_surface_marker_db`

Fallback cache directory:

- `PDOs_outs/Auto_five_state_surface_markers/db_cache`

The script writes a manifest of the actual files used to:

- `PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_database_manifest.csv`

### 3.2 UniProt reviewed human table

The script downloads the current reviewed human UniProtKB table from the REST API with these fields:

- `accession`
- `gene_primary`
- `protein_name`
- `cc_subcellular_location`
- `ft_transmem`
- `ft_topo_dom`

UniProt is the primary source for current subcellular localization and topology.

The script derives these gene-level flags from the UniProt text:

- `uniprot_surface_flag`
  - `TRUE` if subcellular location text contains cell-surface / plasma-membrane phrases such as `cell membrane`, `plasma membrane`, `cell surface`, `apical cell membrane`, `lateral cell membrane`, or `basolateral cell membrane`
- `uniprot_extracellular_flag`
  - `TRUE` if the topological-domain field contains `Extracellular`
- `uniprot_transmem_flag`
  - `TRUE` if the transmembrane feature field contains `TRANSMEM`

### 3.3 ETH Zurich human surfaceome Table S3

The script downloads:

- `https://wlab.ethz.ch/surfaceome/table_S3_surfaceome.xlsx`

The following fields are reused:

- `Surfaceome Label`
- `Surfaceome Label Source`
- `TM domains`
- `topology`
- `CSPA category`
- `UniProt subcellular`

The script derives these gene-level flags:

- `surfaceome_surface_flag`
  - `TRUE` if `Surfaceome Label == "surface"`
- `surfaceome_extracellular_flag`
  - `TRUE` if the surfaceome topology string contains a non-cytoplasmic segment (`NC:`)
- `surfaceome_transmem_flag`
  - `TRUE` if `TM domains > 0`
- `cspa_support_flag`
  - `TRUE` if any `CSPA category` annotation is present
- `cspa_high_confidence_flag`
  - `TRUE` if the CSPA annotation contains `high confidence`

Important: in the implemented workflow, CSPA is treated as supporting evidence only. It is not by itself enough to make a gene a final FACS candidate.

### 3.4 GO cellular-component support

The script uses `org.Hs.eg.db` and `GO.db` locally.

It expands descendants of these GO cellular-component roots:

- `GO:0009986` cell surface
- `GO:0009897` external side of plasma membrane
- `GO:0005886` plasma membrane
- `GO:0005887` integral component of plasma membrane

From these expanded term sets it derives:

- `go_surface_flag`
- `go_plasma_membrane_flag`
- `go_external_side_flag`

These are combined into:

- `go_surface_any_flag = go_surface_flag | go_plasma_membrane_flag`

## 4. Per-State Expression Metrics Used For Sorting

For every retained candidate gene, the script computes state-wise summaries directly from the PDO object.

### 4.1 Target-state metrics

For the assigned target state of each gene:

- `target_pct_cells`
  - median fraction of cells across all valid samples with non-zero counts in that state
- `target_mean_logexpr`
  - median of the mean normalized log-expression across all valid samples in that state

### 4.2 Off-target metrics

Across the other four finalized states, the script identifies:

- `max_off_target_pct_cells`
- `max_off_target_pct_state`
- `max_off_target_mean_logexpr`
- `max_off_target_expr_state`

It then computes:

- `pct_margin = target_pct_cells - max_off_target_pct_cells`
- `expr_margin = target_mean_logexpr - max_off_target_mean_logexpr`

This is the key practical specificity summary for sorting: the gene must not only be present in the target state, but must also beat the strongest off-target state.

## 5. Surface / Topology Gate

### 5.1 Combined annotation flags

The script collapses the annotation sources into:

- `annotation_surface_flag`
  - `TRUE` if the gene has surface/plasma-membrane support from UniProt, GO, or the surfaceome table
- `extracellular_epitope_flag`
  - `TRUE` if UniProt topology, GO external-side annotation, or the surfaceome topology supports extracellular exposure
- `membrane_anchor_flag`
  - `TRUE` if UniProt or surfaceome indicates a transmembrane segment

### 5.2 Final FACS-support definition

The implemented hard gate is:

`supported_for_facs`

This is `TRUE` only if all of the following are true:

- `annotation_surface_flag`
- `extracellular_epitope_flag`
- `membrane_anchor_flag | surfaceome_surface_flag`
- `pct_margin > 0`
- `expr_margin > 0`

This keeps the gate intentionally short:

1. surface/plasma-membrane support
2. extracellular exposure support
3. membrane-anchor or explicit surfaceome support
4. stronger target-state prevalence than any off-target state
5. higher target-state expression than any off-target state

### 5.3 Potential subset

The stricter label is `potential`, which replaces the former `recommended` tier.
This is `TRUE` only if:

- `supported_for_facs == TRUE`
- and either `surfaceome_surface_flag == TRUE`
- or `uniprot_surface_flag == TRUE` together with `membrane_anchor_flag == TRUE`

This gives preference to candidates with stronger membrane / cell-surface support from the external annotation layers.

## 6. Ranking Logic

Ranking is performed within each PDO state.

### 6.1 Continuous terms converted to within-state percent ranks

For each state, the script computes:

- `annotation_rank`
- `target_pct_rank`
- `pct_margin_rank`
- `expr_margin_rank`
- `recurrence_rank`
- `logfc_rank`

These are percent ranks computed within that state's candidate list.

### 6.2 Final ranking score

The current ranking score has been simplified for clarity and rationality. It is:

`facs_priority_score = as.integer(annotation_surface_flag) + target_pct_rank + pct_margin_rank + expr_margin_rank + recurrence_rank`

This explicitly favors:

- Candidates with any valid annotation surface flag (binary +1)
- `target_pct_rank`: broader median target-state cell coverage across samples
- `pct_margin_rank`: better percentage separation from the strongest off-target state
- `expr_margin_rank`: stronger median target-state expression separation
- `recurrence_rank`: higher reproducibility (fraction of eligible samples where the gene was a positive hit)

`annotation_rank` (percentile rank of the discrete `annotation_score`) is no longer included in the priority score to avoid over-complicating the rank, although `annotation_score` itself remains in the table.

### 6.3 Final ordering

Rows are finally ordered by:

1. `state`
2. `recommended_for_facs` descending
3. `supported_for_facs` descending
4. `facs_priority_score` descending
5. `target_pct_cells` descending
6. `pct_margin` descending
7. `expr_margin` descending
8. `sample_recurrence` descending
9. `median_log2FC_hit` descending
10. `gene`

`state_rank` is then assigned from this final order.

## 7. Output Files

The script writes:

- `PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_ranked.csv`
  - full scored table including annotation text, flags, margins, and ranks
- `PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_state_metrics.csv`
  - per-gene per-state `pct_cells` and `mean_logexpr`
- `PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_database_manifest.csv`
  - source URLs, cache paths, and file timestamps for the downloaded references
- `PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_candidates.xlsx`
  - Excel workbook with:
    - `Top5_per_state` (top 5 markers per state with integrated expression/pct data)
    - one sheet per PDO state
- `PDOs_outs/Auto_five_state_surface_markers/Auto_five_state_surface_marker_dotplot.pdf`
  - A bubble plot visualizing the top 5 markers' expression and percent expressed in each state

The Excel workbook is the compact handoff for marker review, while the CSV keeps the full raw annotation text and all scoring columns for auditing or re-ranking.
