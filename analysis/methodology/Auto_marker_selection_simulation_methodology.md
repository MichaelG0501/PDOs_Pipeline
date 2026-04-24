# Auto Marker Selection Simulation

The marker-only classifier does not use finalized state labels during assignment.
All marker-expression summaries and gates use the Seurat RNA `data` layer (log-normalized).
The script validates markers against `Auto_five_state_markers/cache/state_specificity.rds`.

Simulation 1 uses only cells with finalized PDO state labels.
Assignment Logic: For each 3-marker panel, the mean expression of the markers is calculated for each single cell. A panel-specific threshold is set to the 90th percentile of off-target cells or a floor of 0.05 (whichever is higher).
 - If a cell exceeds the threshold for exactly one panel, it is assigned to that state.
 - If a cell exceeds the threshold for multiple panels, it is labeled 'Ambiguous'.
 - If a cell exceeds no thresholds, it is labeled 'Unresolved'.

Simulation 2 creates paired condition shifts (gained/lost states) to test qRT-PCR capturing of relative shifts.
Simulation 3 generates a 10-timepoint longitudinal series to compare state abundance against panel expression (both scaled 0-1) and includes an immune-marker negative control panel.

The immune-regulated panel is used as a negative control; it is not an assignable PDO state in Simulation 1.
