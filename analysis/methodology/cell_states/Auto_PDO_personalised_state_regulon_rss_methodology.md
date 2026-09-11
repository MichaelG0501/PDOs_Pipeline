# Personalised PDO state-regulon RSS methodology

## Purpose and design

This terminal analysis asks whether scATLAS-derived regulons that are strong
and state-specific in a large patient atlas remain active and specific in only
a subset of PDO specimens. Its statistical unit is a PDO specimen (`orig.ident`).
The four patients with matched untreated and FLOT-treated PDOs contribute two
separate specimens because treatment can alter regulon activity; the output
therefore contains 24 specimens representing 20 patients. No patient-level
inference or treatment effect is claimed.

The prespecified tests are:

| Display regulon | Exact scATLAS regulon | Expected PDO state |
| :--- | :--- | :--- |
| SNAI1 | `SNAI1_extended` | Stress-adaptive |
| RXRB | `RXRB` | Stress-adaptive |
| ZBTB14 | `ZBTB14` | Classic proliferation |
| TP73 | `TP73` | Classic proliferation |

`SNAI1_extended` is used because this is the SNAI1 regulon present in the
scATLAS AUCell result. The four-gene non-extended SNAI1 regulon was inferred but
was not present in that scored object. The exact frozen target genes are
embedded in the script and exported to
`tables/Auto_regulon_signature_genes.csv`; this keeps the analysis reproducible
if the large scATLAS SCENIC work cache is removed.

## Cells and state labels

The script intersects cells in `PDOs_merged.rds` with the exactly named current
centred state vector
`centred_mp_refinement/centred_refined_noreg_states.rds`. It uses all cells in
the six canonical states and excludes Hybrid and Unresolved cells from both the
expected-state and reference groups. This prevents ambiguous cells from
weakening or inflating a state contrast. `SUR843T3_PDO` is always excluded in
accordance with repository policy.

This differs deliberately from the existing final-MP SCENIC analysis, whose
5,250 cells were selected and capped by MP assignment for network inference.
Using all canonical-state cells improves specimen-level representation and
avoids conditioning this validation on strong MP assignment.

## Cross-dataset regulon projection

Three requested regulons (SNAI1, RXRB, and ZBTB14) are absent from the
PDO-inferred SCENIC AUC object, while TP73 is present. Their absence from the
PDO network inference is not equivalent to biological inactivity. The workflow
therefore uses all four scATLAS regulons as fixed gene signatures and projects
them into PDO RNA counts with AUCell.

For every canonical-state PDO cell, genes are ranked by expression. Normalized
AUCell AUC is calculated over the top 5% of ranked genes using the target genes
present in the PDO count matrix. A regulon must retain at least five target
genes or the run stops. Per-cell projected scores, target-gene coverage,
metadata, and the AUCell rank parameter are cached in the live output tree so
all tables and figures can be regenerated without rerunning the scoring step.

Absolute AUC values can differ between regulons because target-set size and
composition differ. Absolute activity is therefore compared between specimens
within the same regulon, not between different regulons.

## Per-specimen specificity metrics

SCENIC `calcRSS()` normalizes the regulon's cell-level AUC distribution and
compares it with each state-label distribution by Jensen-Shannon divergence.
For each specimen and regulon, RSS is computed over every canonical state
present in that specimen. The primary specificity gap is:

`next-best RSS gap = RSS(expected state) - max(RSS(each other canonical state))`

This matches the specificity-gap definition in
`Auto_PDO_scAtlas_scenic_comparison.R`. A positive value means the expected
state is the most specific individual state in that specimen; comparison with
the maximum competitor prevents a strong alternative-state association from
being hidden by averaging or pooling weaker states. All six current canonical
PDO states are eligible competitors.

As a sensitivity analysis, the script also combines every non-expected
canonical-state cell into one pooled-rest label and reports the binary RSS gap.
This standardized two-category contrast remains in the source table but does
not drive the primary summary or figures.

The source table additionally reports expected-state and rest mean/median AUC,
their difference, the most active alternative state, and exact cell counts.
These activity summaries are necessary because RSS measures specificity rather
than absolute activity: a positive RSS gap alone does not prove strong activity.

## Reliability and interpretation

The shared repository thresholds of at least 10 expected-state cells and 10
other canonical-state cells are used as a reliability flag. Results below
either count remain in the source tables with their explicit status, but their
figure cells/points are blank so unstable values are not visually promoted.
These thresholds are descriptive safeguards, not formal significance cutoffs.

Cells within one organoid are not independent patient replicates. Accordingly,
the workflow does not perform cell-level hypothesis tests or report p-values.
The analysis is exploratory and personalized: a high activity/high positive-gap
specimen nominates a context in which the atlas-derived regulon is supported,
not a population-level biomarker or treatment recommendation. A low or negative
result can reflect true absence, state scarcity, target dropout, culture
adaptation, or limitations of projecting a reference network across datasets.

## Figures and storage

The evidence matrix shows the primary next-best-state RSS gap on a fixed -1 to
1 color scale; low-cell-support combinations are blank. The multipage profile
figure shows, for each regulon, a dumbbell between expected-state RSS and the
highest-RSS alternative canonical state plus an activity-specificity plot. The
latter identifies specimens that combine positive next-best specificity with
relatively high expected-state AUCell activity without introducing an
arbitrary activity threshold.

All per-cell scores, source tables, figures, and logs are stored under the live
`PDOs_outs/final_mp_scenic/Auto_personalised_state_regulons/` directory.
`PDO_FORCE_REBUILD=1` recomputes AUCell projection, whereas
`PDO_REPLOT_ONLY=1` requires the live cache and regenerates summaries and plots.
The analytical run must use the registered PBS wrapper in the `dmtcp`
environment.
