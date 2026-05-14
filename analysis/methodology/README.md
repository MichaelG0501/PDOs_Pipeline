# PDO Analysis Methodology Index

This folder stores operational methodology notes for `analysis/`. New analysis
scripts must reference a methodology file in their header and update
`analysis/ANALYSIS_MAP.md`.

## Folder Structure

- `shared/` - shared constants, helper functions, output-tier conventions, logging, and cache/replot controls.
- `cell_states/` - state definition, final-state relabeling, marker workflows, SCENIC, FLOT response, drug reversal, and state comparison scripts.
- `clinical/` - clinical association, survival, and presentation figure workflows.
- `metaprograms/` - GeneNMF nMP selection, UCell scoring, metaprogram correlations, and cross-dataset MP analyses.
- `enrichment/` - GO, Hallmark, 3CA, developmental, and pathway enrichment workflows.
- `cnv/` - InferCNA, CNA subclone, Numbat, and CNA diagnostic workflows.
- `demultiplex/` - CellRanger/Souporcell/genotype-assignment rerun workflows.
- `trajectory/` - RNA velocity and pseudotime workflows.
- `plotting/` - general plotting utilities and final-figure style conventions.

## Required Methodology Scope

Each methodology file should document:

1. Scientific aim and current status.
2. Exact input objects, external references, and download requirements.
3. Main algorithmic steps and thresholds.
4. Output files, grouped into `intermediate/`, `tables/`, `figures/`, `logs/`, and `reports/` where applicable.
5. Cache/replot behavior for long-running scripts.
6. Downstream dependencies and scripts that should not use the outputs.
