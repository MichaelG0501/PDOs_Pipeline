# PDO Plotting Methodology

This document describes general plotting conventions for `analysis/plotting/`
and figure-generating scripts across `analysis/`.

## 1. Presentation Requirement

Most PDO downstream figures are used in PowerPoint slides. Plots must remain
readable after insertion into slides. Scripts should therefore choose font
sizes, legend sizes, row/column labels, point sizes, and line widths relative to
the final output dimensions.

## 2. Defaults

New scripts should reuse:

- `analysis/shared/Auto_pdo_analysis_config.R`
- `analysis/shared/Auto_pdo_analysis_helpers.R`

Important defaults:

- canonical state order
- canonical state colors
- slide PDF size `13.333 x 7.5`
- base font size at least 13 for slide figures unless a publication-specific
  figure contract requires smaller text
- PDF primary export
- PNG secondary export only when needed

## 3. Heatmap Rules

For heatmaps:

- avoid unreadably small row names
- rotate column labels only when needed
- size legends relative to the heatmap body
- cache the plotting matrix before drawing if the upstream calculation is slow
- use one figure page per logical panel when labels would otherwise shrink too
  much

## 4. Cache/Replot Rule

Visualization scripts should be able to regenerate figures from cached
`intermediate/` or `tables/` outputs when only aesthetics change. Long-running
data recomputation must not be required for font-size, label, legend, or color
changes.

Recommended environment controls:

- `PDO_REPLOT_ONLY=1`
- `PDO_FORCE_REBUILD=1`

## 5. Legacy Plotting Scripts

Old plotting scripts that use pre-final state vectors, older naming conventions,
or superseded clinical figure layouts should be marked as `legacy` in their
headers and in `analysis/ANALYSIS_MAP.md`.
