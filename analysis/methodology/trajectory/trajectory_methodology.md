# PDO Trajectory And Velocity Methodology

This document describes trajectory scripts under `analysis/trajectory/`.

## 1. Aim

The trajectory workflows study directional relationships between PDO states
using two complementary approaches:

- Monocle3 pseudotime on pre-final four-state noreg state calls
- scVelo RNA velocity using per-sample velocyto outputs and existing PDO UMAPs

## 2. State Inputs

Pseudotime currently uses the pre-unresolved-relabel four-state vector:

- `PDOs_outs/Auto_PDO_states_noreg.rds`

This is deliberate because the trajectory root and state-distance analyses are
defined on the four PDO-intrinsic states:

- `Classic Proliferative`
- `Basal to Intest. Meta`
- `SMG-like Metaplasia`
- `Stress-adaptive`

Velocity metadata exports use both:

- `PDOs_outs/Auto_PDO_states_noreg.rds`
- `PDOs_outs/Auto_PDO_final_states.rds`

## 3. Pseudotime Script Roles

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

## 4. Velocity Script Roles

The velocity workflow exports PDO metadata and QC barcodes, prepares references,
runs per-sample CellRanger/velocyto/scVelo, and builds state-direction node
plots. Most velocity scripts are currently untracked and should remain unstaged
unless explicitly requested.

## 5. Output Standards

Trajectory scripts should cache heavy objects in `intermediate/` and write
final state-distance tables and plots to `tables/` and `figures/`. Replot-only
mode should regenerate PDFs from cached per-sample trajectory assets.
