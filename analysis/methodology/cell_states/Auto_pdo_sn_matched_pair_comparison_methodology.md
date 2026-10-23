# Matched PDO and snRNA-seq comparison using current centred inputs

## Scope

The retained entry-point basename is
`analysis/cell_states/legacy_Auto_pdo_sn_matched_pair_comparison.R`; its current
registry marks it active terminal. This user-authorized update replaces its
historical analytical route. Compare SUR680T3_PDO with H_post_T1_biopsy and
SUR791T3_PDO with L_post_T1_biopsy. Exclude SUR843T3_PDO as in other PDO analyses.
The comparison is descriptive for two matched patients; no cell-level hypothesis
tests or population-level inferences are made.

## Canonical inputs

Read live PDOs_outs/PDOs_merged.rds and the current PDO centred_mp_refinement
MP genes, UCell scores, named noreg states, adjusted score matrix, and grouping
CSV, exactly as documented in PDO AGENTS.md. Read live snSeq_Pipeline/sn_outs
malignant metadata and the corresponding files under cell_states/centred/
intermediate and tables. snRNA scores use fixed current centred scRef MPs; native
PDO MPs were independently derived. Never use uncentred nMP13/nMP19 objects,
historical Auto_final_states objects, 3CA state relabelling, or path fallbacks.

## Denominators and MP interpretation

Use exact metadata rownames/barcodes and orig.ident sample identities. Require
every matched cell to have a known state, all grouped adjusted MP scores, and all
UCell scores. Missing matched samples, cells, MPs, or non-finite scores stop the
run. Preserve native state definitions: PDO has six biological states; projected
scRef/snRNA has five. Do not rename biologically related labels into identical
states or treat the same MP ID in the two collections as the same programme.

Use the producer's persistent adjusted score matrices directly, without
recentring a matched subset. The top-MP table contains both the highest adjusted
MP including CC/QC rows and the highest adjusted non-CC state-defining MP.
Deterministic ties follow grouping-table order. All matched cells, including
Hybrid and Unresolved cells, remain in both state and MP denominators. Complete
missing state/MP categories with zero counts. Each sample's percentages sum to
100 within an analysis. No historical overlap/correlation metric is retained:
the native state spaces and independent MP identifiers differ.

## Outputs and execution

Write a per-cell assignment table, zero-completed state/MP proportion tables,
sample-level cell counts and dominant states, a persistent result cache, figure
index, run log, and one two-page PDF (one patient per page). State and non-CC
MP panels remain separate for each modality. All production outputs live under
PDOs_outs/Auto_pdo_sn_matched_pair_comparison/{intermediate,tables,figures,logs,reports}.
No stale output is deleted automatically.

Submit job_scripts/cell_states/legacy_Auto_pdo_sn_matched_pair_comparison.sh:
dmtcp, 2 CPUs, 32 GB, 2 hours. Default rebuilds cached source tables; replot with
PDO_REPLOT_ONLY=1 only if input paths, sizes, and modification times match. This
fingerprint avoids hashing the 2.5-GB merged object on every replot; it cannot
detect content changes preserving both size and mtime. SNSEQ_OUTS_DIR and
PDO_SN_COMPARISON_OUT_DIR allow explicitly supplied isolated validation roots;
production defaults are the live canonical paths. Smoke-test outputs are
synthetic validation artifacts and must never be presented as biological results.
