# Centred MP Ordering, State Definition, and Downstream Scoring

## Scope and status

This is the active methodology for centred steps 05–08:

- `Auto_05_centred_refined_mp_ordered_heatmaps.R`
- `Auto_06_centred_refined_state_definition_noreg.R`
- `Auto_07a_ucell_scoring.R`
- `Auto_07_3ca_vs_refined_mp_correlation.R`
- `Auto_08_tcga_mp_survival_volcano_centred.R`

The scripts consume the filtered live outputs of step 04. They do not restore
MPs removed by parent/refined-MP QC, and they do not use the superseded
uncentred nMP=13 objects or `Auto_PDO_final_states.rds` route.

## Authoritative MP set and ordering

Step 04 is the only final MP filtering point. It retains a refined MP when it
has at least five genes and NMF-program coverage in at least three observed PDO
samples. The live `merged_refined_mp_genes.rds` is authoritative. Step 05
intersects the assignment table with this gene-list object and applies no
additional coverage or gene-count cutoff.

The current biological display order is:

1. cell cycle: MP11, MP1, MP2, MP3;
2. Classic proliferation: MP19+;
3. Basal to intestinal metaplasia: MP15, MP5+, MP12;
4. SMG to intestinal metaplasia: MP13b, MP14b, MP16b, MP17+, MP8+;
5. Stress adaptive: MP9;
6. PDO medium induced: MP18.

Step 05 writes this mapping to
`tables/centred_refined_mp_state_grouping.csv` and the complete order, colors,
descriptions, and group mapping to `centred_refined_mp_strict_order.rds`. These
are persistent live inputs for later plotting and matched-treatment workflows.

Program-similarity order is inherited from the optimal centred GeneNMF tree
within each final MP. Score correlation uses Spearman correlation. The split
correlation panel preserves the sub-MP identities from step 03 to audit which
components were merged by step 04.

## UCell scoring

Step 07a scores the final PDO MP gene lists and the external 3CA signatures on
the log-normalized RNA expression matrix with `ScoreSignatures_UCell`. UCell is
rank based; the script removes the generated `_UCell` suffix and stores cells
as rows and signatures as columns. Final PDO scores are written to both live
and ephemeral storage. The live copies are authoritative downstream inputs.

## Noreg state definition

Step 06 uses only non-cell-cycle groups for state assignment. For each MP, it:

1. subtracts the within-sample mean UCell score;
2. divides by that MP's standard deviation in the corresponding PDO batch;
3. replaces non-finite adjusted values with zero;
4. takes the maximum adjusted MP score within each biological state group.

For a cell with group maxima \(g_1 \ge g_2\):

- `Unresolved` when \(g_1 < 0.5\);
- `Hybrid` when \(g_1 \ge 0.5\) but \(g_1-g_2 < 0.3\);
- otherwise, the state belonging to \(g_1\).

The 0.5 activity threshold and 0.3 separation threshold are fixed values
inherited from the noreg Approach-B analogue. They are not optimized on the
current outcome data and are not varied by sample. Cell-cycle MPs remain in QC
heatmaps but cannot define a state.

The state vector, adjusted MP matrix, and group-max matrix retain exact cell
barcode names and are saved to both live and ephemeral storage. Those live RDS
files are required to reproduce state figures without rerunning UCell or NMF.

## 3CA comparison

Step 07 compares cell-level 3CA and current PDO centred-refined scores using
Spearman correlation. Its unresolved-cell panels are descriptive: they show
which external programs are active among cells that fail the state threshold,
but do not relabel those cells and do not create a downstream state vector.

## TCGA survival analysis

Step 08 reads persistent current MP gene lists, TCGA ESCA TPM, and reconstructed
clinical metadata. It filters signatures to genes present in the bulk matrix,
requires at least five genes per signature, computes sample-level gene-set
scores using the GSVA implementation in the script, and reports Cox-model and
split-method sensitivity results. These are association analyses, not causal
or treatment-response estimates. The CSV contains the tested variants; the PDF
is a terminal visualization.

## Run order and environments

After step 04 completes successfully:

1. submit `Auto_run_centred_07a_ucell.sh` in `gnmf` so the final PDO and 3CA
   score matrices are refreshed from exactly the step-04 gene lists;
2. submit `Auto_replot_05.sh` in `dmtcp` after scoring;
3. submit `Auto_centred_06.sh` after step 05 so state calls use that refreshed
   score matrix and the persisted step-05 grouping/order;
4. submit `Auto_run_centred_07_08_terminal.sh` after step 06 for the 3CA
   comparison and TCGA survival outputs.

All analytical execution must use PBS with `#PBS -koed`. The scripts prefer
live inputs and use ephemeral copies only as backward-compatible caches.
