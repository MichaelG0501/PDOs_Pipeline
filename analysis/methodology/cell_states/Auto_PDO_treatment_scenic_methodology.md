# Matched PDO Treatment SCENIC Methodology

## Design

`Auto_PDO_treatment_scenic.R` analyzes exactly four patient-matched PDO pairs:
SUR1070, SUR1072, SUR1090, and SUR1181, each with untreated and FLOT-treated
samples. All eight required samples must be present. The script does not use a
cell-state vector, so centred-state changes do not alter its regulon inference.

## SCENIC input and filtering

Raw RNA counts from `PDOs_merged.rds` are restricted to the eight samples.
Genes pass SCENIC filtering when they satisfy the script's adaptive minimums:

- total counts at least `max(3 * 0.01 * number_of_cells, 20)`;
- detected samples/cells at least `max(0.01 * number_of_cells, 20)`.

These bounds follow SCENIC's count/detection logic while enforcing an absolute
minimum of 20 so extremely sparse genes do not enter GENIE3. Human hg38
cisTarget ranking databases are supplied explicitly through `db_dir`. The
workflow runs co-expression inference, motif pruning, AUCell regulon scoring,
and regulon-specificity scoring. SCENIC databases and annotation patches are
version-sensitive; the run log and environment must be retained with results.

## Treatment comparisons

Cell-level regulon AUC is summarized by treatment and by sample. The main
treatment screen uses a two-sided Wilcoxon test on cell-level AUC and applies
Benjamini-Hochberg correction across regulons; the default display threshold is
FDR 0.05 and can be changed only through the documented `fdr_threshold`
argument.

Because cells from one PDO are not independent patients, the workflow also
computes treated-minus-untreated mean AUC separately for each of the four
patients. Directional consistency records how many patient pairs increase or
decrease. With only four pairs, these paired changes are descriptive and should
take precedence over isolated cell-level significance when interpreting
treatment effects.

For network figures, regulon targets are retained above the within-regulon
median importance score. This 50th-percentile cutoff is a visualization rule,
not a significance threshold or evidence that lower-ranked targets are absent.

## Storage and reruns

Final tables, figures, AUC matrices, regulons, and RSS matrices are written
under live `PDOs_outs/treatment_scenic/`. Heavy SCENIC work/cache files
are kept under the matching ephemeral path. `prepare_only=true` validates the
sample design without heavy computation. A full run must use the PBS wrapper
with `#PBS -koed`; it must never run on a login node.
