# PDO Trajectory Methodology

This document describes the trajectory analysis scripts under `analysis/trajectory/`.
The two complementary approaches — **Monocle3 pseudotime** and **scVelo RNA
velocity** — are documented in separate sections.

---

## 1. Common State Inputs

Both workflows operate on the four PDO-intrinsic malignant states:

| State | Color |
|---|---|
| Classic Proliferative | `#E41A1C` |
| Basal to Intest. Meta | `#4DAF4A` |
| SMG-like Metaplasia | `#FF7F00` |
| Stress-adaptive | `#984EA3` |

Source files:

- **Pre-unresolved-relabel four-state vector** (`state_four`):
  `PDOs_outs/Auto_PDO_states_noreg.rds`
- **Finalized five-state vector** (`state_final`):
  `PDOs_outs/Auto_PDO_final_states.rds`

Pseudotime uses `state_four` exclusively so that root placement and distance
calculations are defined on the four PDO-intrinsic states.

Velocity metadata exports both vectors: `state_four` for transition-direction
summaries, `state_final` for UMAP coloring.

---

# Part A — Pseudotime (Monocle3)

## A1. Script Roles

### `Auto_PDO_pseudotime_helpers.R`

Shared helper library for building per-sample Monocle3 trajectories and state
distance summaries.

### `Auto_PDO_pseudotime_samples.R`

Builds and caches per-sample trajectory objects, pseudotime vectors, and
metadata tables.

### `Auto_PDO_pseudotime_linear_plot.R`

Reuses cached trajectory assets to regenerate presentation reports without
recomputing trajectories.

### `Auto_PDO_pseudotime_state_distance_matrix.R`

Computes state-to-state distances from cached pseudotime, principal graph, and
UMAP coordinates.

## A2. Output Standards

Trajectory scripts cache heavy objects in `intermediate/` and write
final state-distance tables and plots to `tables/` and `figures/`. Replot-only
mode regenerates PDFs from cached per-sample trajectory assets.

## A3. Pseudotime PDF Outputs

| PDF File | Script |
|---|---|
| `Auto_PDO_pseudotime_combined.pdf` | `Auto_PDO_pseudotime_samples.R` |
| `Auto_PDO_pseudotime_linear_reports.pdf` | `Auto_PDO_pseudotime_linear_plot.R` |
| `Auto_PDO_state_distance_method_comparison_heatmap.pdf` | `Auto_PDO_pseudotime_state_distance_matrix.R` |
| `Auto_PDO_state_distance_nodeplot.pdf` | `Auto_PDO_pseudotime_state_distance_matrix.R` |

---

# Part B — RNA Velocity (scVelo)

## B1. Pipeline Overview

The RNA velocity pipeline consists of six scripts executed in a strict dependency
chain. The master submission script (`Auto_00_submit_pdo_velocity.sh`)
orchestrates the entire sequence via PBS `qsub` with job-dependency chaining.

### Execution Order

```
Step 0  Auto_export_pdo_velocity_metadata.R    (dmtcp)   — export cell metadata, barcode lists, manifest
Step 0  Auto_prepare_pdo_velocity_refs.py       (velocity) — decompress gene GTF, download + convert RepeatMasker GTF
Step 1  Auto_run_pdo_cellranger.sh              (per-sample, Cynthia batch only) — re-run CellRanger count to produce BAMs
Step 2  Auto_filter_sort_pdo_velocity.sh        (per-sample) — filter BAM to QC barcodes, coordinate-sort
Step 3  Auto_run_pdo_velocyto.sh                (per-sample) — run velocyto to produce per-sample loom files
Step 4  Auto_scvelo_pdo_visualise.py            (all samples) — run scVelo, produce all PDF visualisations
Step 5  Auto_pdo_velocity_nodeplots.R           (all samples) — R-based treatment-comparison node plots
```

### Master Orchestrator — `Auto_00_submit_pdo_velocity.sh`

Runs Steps 0 interactively, then loops over the sample manifest and submits
Steps 1-3 as PBS jobs with `afterok` dependencies. After all velocyto jobs
complete, submits Step 4 as a single dependent job. Step 5 is run independently.
Throttle limit: 46 concurrent PBS jobs.

---

## B2. Step 0 — Metadata Export

### Script: `Auto_export_pdo_velocity_metadata.R`

**Inputs**: `PDOs_merged.rds`, `Auto_PDO_states_noreg.rds`,
`Auto_PDO_final_states.rds`.

**Outputs** (all under `PDOs_outs/Auto_velocity_PDO/`):

| File | Content |
|---|---|
| `tables/Auto_pdo_velocity_cell_metadata.csv` | Per-cell metadata: cell_id, sample, raw_barcode, batch_type, treatment, state_four, state_final, umap_1, umap_2, nCount_RNA, nFeature_RNA |
| `tables/Auto_pdo_velocity_sample_manifest.csv` | Per-sample manifest: sample, batch_type, treatment, pool, fastq_dir, cellranger_out, bam, barcodes_file, n_cells, has_bam |
| `tables/Auto_pdo_velocity_state_four_summary.csv` | Cross-tabulation of cells per sample × state_four |
| `barcodes/<sample>_qc_barcodes.tsv` | Sorted unique raw barcodes that passed QC, one per line |

**Processing details**:

- State vectors are coerced to character with `names()` explicitly restored to
  prevent accidental name loss.
- `raw_barcode` is extracted by stripping the `<sample>_` prefix from
  `cell_id`. If the cell_id does not start with the sample prefix, it falls
  back to stripping everything up to the last underscore.
- `batch_type` is classified as `New_batch` when the sample name matches
  `_(Untreated|Treated)_PDO$`, otherwise `Cynthia_batch`.
- `treatment` is `Untreated`, `Treated`, or `Cynthia` (legacy samples with no
  paired treatment design).
- UMAP coordinates are taken from the `umap` reduction in `PDOs_merged.rds`.

---

## B3. Step 0 — Reference Preparation

### Script: `Auto_prepare_pdo_velocity_refs.py`

Produces two reference files required by velocyto:

1. **Gene GTF**: decompresses the CellRanger GRCh38-2024-A `genes.gtf.gz` into
   `ref/genes.GRCh38-2024-A.gtf`.
2. **RepeatMasker GTF**: downloads UCSC hg38 `rmsk.txt.gz`, converts each
   repeat element to a single-exon GTF line with `gene_id`, `transcript_id`,
   `rep_class`, and `rep_family` attributes. Saved as
   `ref/repeatmasker.hg38.gtf`.

Both files are idempotent (skipped if the output already exists and is
non-empty).

---

## B4. Step 1 — CellRanger Re-run (Cynthia Batch Only)

### Script: `Auto_run_pdo_cellranger.sh`

Only Cynthia-batch samples lack pre-existing BAMs. For these samples the script
runs `cellranger count` (v9.0.1) with `--create-bam=true` against the
GRCh38-2024-A transcriptome to produce `possorted_genome_bam.bam`. New-batch
samples already have BAMs from the demultiplexing pipeline and this step exits
immediately.

**PBS resources**: 16 cores, 512 GB RAM, 48 h walltime.

---

## B5. Step 2 — BAM Filtering and Coordinate Sorting

### Script: `Auto_filter_sort_pdo_velocity.sh`

For each sample:

1. Reads the input BAM (CellRanger output).
2. Filters alignments: only reads whose `CB:Z:` tag matches a barcode in the
   QC barcode list are retained. Header lines are passed through.
3. Pipes directly into `samtools sort` for coordinate-order output.
4. Indexes the output BAM and runs `samtools quickcheck`.

**Output**: `coord/<sample>.qc.coord.bam` + `.bai` index.

**PBS resources**: 8 cores, 96 GB RAM, 24 h walltime.

---

## B6. Step 3 — Velocyto Loom Generation

### Script: `Auto_run_pdo_velocyto.sh` → `Auto_velocyto_pdo_run.py`

The shell wrapper activates the `velocity` conda environment and calls the
Python runner.

### `Auto_velocyto_pdo_run.py` — Monkey-patching velocyto

The velocyto library's internal `ExInCounter.peek()` method is replaced with a
custom version that correctly detects CellRanger BAM tags (`CB` + `UB`). The
original velocyto code can fail when these tags are present but the BAM is
coordinate-sorted (velocyto was designed for 10x position-sorted BAMs with the
older `CR` tag convention).

The monkey-patched `_peek_cellranger()` inspects the first 1000 mapped reads
and sets:
- `cellbarcode_str = "CB"`, `umibarcode_str = "UB"` (CellRanger)
- or `"CB"` + `"pN"` (Parse Biosciences)
- or `"XC"` + `"XM"` (Drop-seq)

The patched `_peek_umi_only()` does the same for UMI-only detection.

After patching, `velocyto.commands._run._run()` is called with:
- `logic="Default"` (standard exon/intron counting)
- `multimap=False`
- `without_umi=False`
- RepeatMasker GTF as the mask file (to exclude repetitive elements)

**Output**: `looms/<sample>/<sample>.loom` — a loom file containing spliced,
unspliced, and ambiguous count layers.

**PBS resources**: 8 cores, 96 GB RAM, 36 h walltime.

---

## B7. Step 4 — scVelo Analysis and Per-Sample Visualisation

### Script: `Auto_scvelo_pdo_visualise.py`

This is the core velocity analysis and visualisation script. It loads each
sample's loom file, runs scVelo, and produces three multi-page PDFs plus CSV
summary tables.

**PBS resources**: 8 cores, 160 GB RAM, 24 h walltime.

### B7.1 Sample Loading (`load_sample`)

For each sample:

1. Reads the single `.loom` file from `looms/<sample>/`.
2. Extracts raw barcodes from loom observation names by splitting on `:` and
   taking the last token, then stripping any trailing `x`.
3. Matches loom barcodes to the QC barcode list from the exported metadata CSV.
   If a barcode is missing and doesn't end with `-1`, the `-1` suffix is
   appended and retried. This handles the CellRanger convention where 10x
   barcodes carry a `-1` suffix.
4. Subsets the AnnData to only QC-passing cells and renames observations to
   Seurat `cell_id` identifiers.
5. Transfers all metadata columns (sample, batch_type, treatment, state_four,
   state_final, etc.) into `adata.obs`.
6. Fills missing state labels with `"Unassigned"`.
7. Inserts the original Seurat UMAP coordinates into `adata.obsm["X_umap"]`.

### B7.2 Velocity Computation (`run_velocity`)

Per sample, the following scVelo + scanpy pipeline runs:

| Step | Function | Parameters |
|---|---|---|
| Gene filter | `scv.pp.filter_genes` | `min_shared_counts=20` |
| Normalise | `scv.pp.normalize_per_cell` | library-size normalisation |
| Log-transform | `sc.pp.log1p` | natural log(1+x) |
| HVGs | `sc.pp.highly_variable_genes` | top 3000, Seurat flavour; skipped if < 3000 genes remain |
| PCA | `sc.tl.pca` | `svd_solver="arpack"`, `n_comps` = min(30, n_obs−1, n_vars−1, ≥2) |
| Neighbours | `sc.pp.neighbors` | `n_neighbors` = min(30, n_obs−1, ≥5), using computed PCs |
| Moments | `scv.pp.moments` | first and second moments of spliced/unspliced counts |
| Velocity | `scv.tl.velocity` | `mode="stochastic"` — stochastic model of transcriptional dynamics |
| Velocity graph | `scv.tl.velocity_graph` | cell-cell transition probability matrix |
| Velocity embedding | `scv.tl.velocity_embedding` | `basis="umap"` — projects velocity vectors into the Seurat UMAP space |
| Confidence | `scv.tl.velocity_confidence` | per-cell velocity length and confidence scores |

**Caching**: results are saved as `h5ad/Auto_scvelo_<sample>.h5ad`. On re-run,
if the h5ad exists and contains the required fields, the cached object is loaded
instead of recomputing.

### B7.3 Sanitisation (`sanitise_velocity_fields`)

After velocity computation (or cache load), all NaN/Inf values in `X_umap`,
`velocity_umap`, `velocity_length`, and `velocity_confidence` are replaced with
0.0 to prevent matplotlib rendering failures.

### B7.4 Display Basis Compression (`prepare_display_basis`)

To create visually comparable per-sample UMAP panels:

1. **Center**: subtract the per-sample median UMAP coordinate.
2. **Compress outliers**: compute the 90th percentile radial distance (`radius`).
   Points within the radius keep their original distance; points beyond are
   compressed using `radius + 0.15 × (distance − radius)`.
3. **Scale**: normalise so that the 99th percentile span maps to ±6 coordinate
   units.
4. Store as `X_umap_display` and `velocity_umap_display` (velocity vectors are
   scaled by the same factor).

### B7.5 Core PDF — Per-Sample Visualisations

**Output**: `figures/Auto_pdo_velocity_per_sample_visualisations.pdf`

Each sample gets one page with **three panels** (left to right):

| Panel | Content | Coloring | Method |
|---|---|---|---|
| 1 — Finalized states | Scatter UMAP | `state_final` (all 5+ states) | `sc.pl.embedding` on `X_umap_display` |
| 2 — Four-state call | Scatter UMAP | `state_four` (4 major + Unresolved/Hybrid/Unassigned) | `sc.pl.embedding` on `X_umap_display`, legend hidden |
| 3 — Velocity stream | Streamlines overlaid on scatter | `state_final` colors, legend hidden | `scv.pl.velocity_embedding_stream` on `X_umap_display` |

- Point size scales inversely with cell count: `max(10, min(40, 26000/n_obs))`.
- Stream-plot points are 80% of the scatter size (minimum 8).
- Pages are rasterised at 190 DPI to control PDF file size.

### B7.6 Extended PDF — Per-Sample QC Metrics

**Output**: `figures/Auto_pdo_velocity_per_sample_visualisations_extended.pdf`

Each sample gets one page with **three panels**:

| Panel | Content | Coloring | Method |
|---|---|---|---|
| 1 — Velocity grid | Arrow grid overlaid on scatter | `state_final` | `scv.pl.velocity_embedding_grid` with `arrow_size=3.8`, `arrow_length=6.5`, `density=0.65`, `scale=0.35`, black arrows |
| 2 — Velocity length | Scatter UMAP | Continuous viridis colormap | `velocity_length` — magnitude of the velocity vector per cell |
| 3 — Velocity confidence | Scatter UMAP | Continuous viridis colormap | `velocity_confidence` — correlation of the cell's velocity with the velocities of its neighbours |

### B7.7 State-Direction Calculation (`state_direction_tables`)

This is the core quantitative summary of state transitions. For each sample:

#### Node Metrics

For each of the four major states:
- **cells**: number of cells assigned to that state (from `state_four`)
- **pct_major**: percentage of cells relative to total cells in the four major
  states only (excludes Unresolved/Hybrid/Unassigned)

#### Edge Metrics — Velocity Alignment Score

For every ordered pair of states (source → target), a **velocity alignment
score** is computed:

1. **Mean UMAP velocity** of the source state: average of `velocity_umap`
   vectors across all cells in the source state (requires ≥ 5 cells).
2. **Target direction**: vector from the source state's UMAP centroid to the
   target state's UMAP centroid.
3. **Alignment** = cosine similarity between the mean velocity vector and the
   source-to-target direction vector:

   ```
   alignment = dot(mean_velocity, target_direction) / (‖mean_velocity‖ × ‖target_direction‖)
   ```

   Range: [−1, +1]. Positive values indicate the source state's cells are
   moving toward the target state in UMAP space. Negative values indicate
   movement away from the target.

Additional edge columns: `source_velocity_norm`, `source_to_target_distance`,
`source_cells`, `target_cells`.

### B7.8 Direction PDF — State-Transition Network Diagrams

**Output**: `figures/Auto_pdo_velocity_state_directions.pdf`

This multi-page PDF contains directed network diagrams showing velocity-based
state transitions.

#### Page Structure

| Page | Content |
|---|---|
| 1 | **All PDOs** (20 samples aggregated) |
| 2 | **Batch comparison** — three panels: Cynthia batch, New batch untreated, New batch treated |
| 3+ | **Paired patients** — one page per patient with both Untreated and Treated samples, side by side |
| Remaining | **Unpaired samples** — one panel per sample |

#### Aggregation Method (`aggregate_direction`)

When multiple samples are combined (e.g. all PDOs, or a batch):

- **Node pct_major**: mean of per-sample percentages
- **Edge velocity_alignment**: mean of per-sample alignment scores
- Additional edge columns: `median_alignment`, `positive_fraction` (fraction of
  samples where alignment > 0), `sample_n` (number of contributing samples)

#### Network Diagram Rendering (`draw_direction_network`)

- **Node positions**: fixed layout — Classic Proliferative top-left (−1, 0.72),
  Basal to Intest. Meta top-right (1, 0.72), SMG-like Metaplasia bottom-right
  (1, −0.72), Stress-adaptive bottom-left (−1, −0.72).
- **Node size**: scales linearly with `pct_major`, range 320–1800 scatter units.
- **Node color**: canonical state colors.
- **Node label**: abbreviated state name + percentage (e.g. "Classic\nProlif.\n32.5%"),
  rendered in bold with white stroke outline.
- **Edge threshold**: only edges with `velocity_alignment ≥ 0.35` are drawn.
- **Edge rendering**: `FancyArrowPatch` with `arc3` curvature (rad = ±0.18,
  direction depends on state ordering in the canonical list).
- **Edge width**: `0.8 + 4.8 × alignment` — thicker arrows for stronger
  transitions.
- **Edge alpha**: `0.18 + 0.65 × alignment` — more opaque for stronger
  transitions.
- **Edge label**: alignment score displayed at midpoint as "0.XX" in a white
  background box.
- If no edge exceeds the threshold, a "No state-pair alignment ≥ 0.35" message
  is displayed.

### B7.9 CSV Table Outputs

All saved to `PDOs_outs/Auto_velocity_PDO/tables/`:

| File | Content |
|---|---|
| `Auto_pdo_velocity_scvelo_cell_metadata.csv` | Per-cell metadata + velocity_umap_1/2 for all samples |
| `Auto_pdo_velocity_state_nodes.csv` | Per-sample, per-state: cells, pct_major |
| `Auto_pdo_velocity_state_direction_edges.csv` | Per-sample, per-source-target: velocity_alignment + auxiliary metrics |
| `Auto_pdo_velocity_top_state_direction_per_source.csv` | Highest-alignment target for each source state in each sample |
| `Auto_pdo_velocity_direction_audit_summary.csv` | Per group (All/Cynthia/New untreated/New treated), per target: top_positive_count, mean/median alignment |
| `Auto_pdo_velocity_group_state_direction_edges.csv` | Aggregated edges for All PDOs + batch groups |

---

## B8. Step 5 — Treatment-Comparison Node Plots (R)

### Script: `Auto_pdo_velocity_nodeplots.R`

**Output**: `figures/Auto_pdo_velocity_nodeplot_untreated_vs_treated.pdf`

This R script reads the CSV tables produced by Step 4 and creates ggplot2
network diagrams comparing untreated vs treated PDOs.

### B8.1 Data Preparation

- Reads `Auto_pdo_velocity_state_nodes.csv` and
  `Auto_pdo_velocity_state_direction_edges.csv`.
- Extracts patient ID from sample name: `sub("^(SUR[0-9]+).*", "\\1", sample)`.
- Uses the same fixed node layout as the Python script.
- Plotting threshold for edges: `velocity_alignment > 0.10` (lower than the
  Python script's 0.35 threshold).

### B8.2 Aggregation

- **Nodes** (`combine_nodes`): median of `cells` and `pct_major` across samples.
- **Edges** (`combine_edges`): mean of `velocity_alignment` across samples,
  then filtered to alignment > 0.10.

### B8.3 Rendering (`draw_nodeplot`)

Each panel is a `ggplot` with:

- `geom_curve` for edges — curvature 0.18, closed arrows, linewidth and alpha
  both scaled by velocity_alignment. Linewidth range: 0.5–5.8. Alpha range:
  0.25–0.9.
- `geom_label` at edge midpoints showing the alignment score to 2 decimal places.
- `geom_point` for nodes — size scaled by pct_major (range 9–30), colored by
  canonical state colors.
- `geom_text` for node labels — abbreviated state name + percentage, bold.
- Shared scale limits across panels (max_node, max_edge) so Untreated and
  Treated panels are visually comparable.

### B8.4 PDF Page Structure

| Page | Content |
|---|---|
| 1 | **New Batch overview** — Untreated (median nodes, mean arrows) vs Treated, side-by-side via patchwork |
| 2+ | **Per-patient pairs** — one page per patient with both Untreated and Treated samples |
| Remaining | **Cynthia batch singletons** — one panel per Cynthia sample |

---

## B9. Velocity Constants and Thresholds Summary

| Parameter | Value | Where Used |
|---|---|---|
| Velocity mode | `stochastic` | `scv.tl.velocity` |
| Min shared counts | 20 | `scv.pp.filter_genes` |
| Top HVGs | 3000 | `sc.pp.highly_variable_genes` |
| PCA components | min(30, n_obs−1, n_vars−1), ≥2 | `sc.tl.pca` |
| Neighbours | min(30, n_obs−1), ≥5 | `sc.pp.neighbors` |
| UMAP outlier compression | 90th percentile radius, 15% compression beyond | `prepare_display_basis` |
| Direction edge threshold (Python) | 0.35 | `draw_direction_network` |
| Direction edge threshold (R) | 0.10 | `Auto_pdo_velocity_nodeplots.R` |
| Min cells for state centroid | 5 | `state_direction_tables` |
| Raster DPI | 190 | `save_raster_page` |
| Point size formula | max(10, min(40, 26000/n_obs)) | `point_size` |
