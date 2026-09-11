# Centred PDO Metaprogram Refinement Methodology

## Status and scope

This is the active methodology for
`analysis/metaprograms/centred/Auto_02_nmf_rank_selection_diagnostics.R`,
`Auto_03_mp_refinement_submp.R`, and
`Auto_04_mp_refinement_merge_correlated_submps.R`. It adapts the current
malignant scRef centred workflow to purely malignant PDO data. No cell-type or
CNA filtering is performed, and `SUR843T3_PDO` remains excluded.

## Run order

1. `Auto_01_centred_geneNMF.R` creates centred multiNMF programs and nMP 4:25
   metaprogram objects.
2. `Auto_02_nmf_rank_selection_diagnostics.R` calculates mean silhouette and
   cosine-distance WSS, selects the silhouette kneedle optimum, runs initial
   enrichment, and plots the optimal program-similarity matrix.
3. `Auto_03_mp_refinement_submp.R` applies parent-MP QC and subdivides MPs with
   silhouette between 0 and 0.2.
4. `Auto_04_mp_refinement_merge_correlated_submps.R` merges correlated sibling
   sub-MPs, recalculates coherence metrics, applies final refined-MP QC, and
   writes the persistent final gene lists, UCell scores, tables, and figures.

The PBS wrappers `Auto_centred_02.sh`, `Auto_centred_03.sh`, and
`Auto_centred_04.sh` use the `dmtcp` environment and must be submitted with
`afterok` dependencies in this order. Set `PDO_FORCE_REBUILD=TRUE` for a full
refinement regeneration.

## Parent-MP QC before splitting

The number of samples is inferred from NMF program names rather than hard-coded.
GeneNMF fractional `sampleCoverage` is converted to an integer sample count.
The current scRef boundary rules are retained literally:

- keep: silhouette at least 0.2, coverage in at least 3 samples, and more than
  5 genes;
- split: silhouette above 0 and below 0.2, coverage in at least 3 samples, and
  more than 5 genes;
- remove: negative silhouette, coverage below 3 samples, or at most 5 genes.

At the current PDO nMP 20 solution this removes MP10 (4 genes) and MP20
(2-sample coverage), in addition to MPs with negative silhouette, before any
sub-MP splitting.

## Splitting and merging

Intermediate-silhouette MPs are hierarchically split using the NMF-program
cosine-similarity matrix. Sub-MP gene signatures use the GeneNMF-style weighted
consensus implementation in step 03. Step 04 merges qualifying siblings when
their mean Spearman correlation exceeds `PDO_SUBMP_MERGE_COR` (default 0.4)
with at least `PDO_SUBMP_MERGE_FRACTION` (default 0.25) of their siblings.

## Final refined-MP QC

After correlated merging, sample coverage is recomputed directly from the
program-to-refined-MP assignments. A refined MP is retained only when it has:

- coverage in at least 3 samples; and
- at least 5 genes.

There is no PDO-specific manual MP removal. The filtering audit is written to
`PDOs_outs/centred_mp_refinement/tables/merged_refined_mp_final_filtering.csv`.
Filtered gene lists, weights, and UCell scores are written to both the
ephemeral cache and live persistent storage. The live copies are authoritative
downstream inputs.

Step-03 inputs needed by step 04 (`split_results.rds`, refined genes, weights,
assignments, and refined UCell scores) are likewise saved at the live centred
output root as well as in the ephemeral cache. Step 04 prefers those live files
and uses the ephemeral paths only as a compatibility fallback for older runs.

## Customized optimal-nMP heatmap

Step 02 follows the current scRef layout but annotates `batch` instead of
`study`. Program sample names are mapped as follows:

- samples ending in `_Treated_PDO` or `_Untreated_PDO`: `batch2`;
- the Cynthia cohort (other historical `SUR...PDO` samples): `pdo`;
- `TEMP_new4samples_*` and finalized SUR1346/SUR1363/SUR1384/SUR1391 names:
  `new4samples`.

The figure is
`figures/Auto_centred_nMP_{optimal}_custom_heatmap_by_batch.pdf`; its exact
program/sample/batch annotations are saved in
`tables/Auto_centred_nMP_{optimal}_program_batch_annotations.csv`.

## Persistent outputs and downstream use

The primary downstream inputs under `PDOs_outs/centred_mp_refinement/` are
`optimal_nMP.rds`, `merged_refined_mp_genes.rds`,
`merged_refined_mp_gene_weights.rds`, and
`merged_refined_ucell_scores.rds`. These live files must be sufficient for
ordered heatmaps, enrichment, state definition, and replotting if ephemeral
storage is removed. Figures produced by step 02 are terminal QC; filtered step
04 objects are active upstream inputs.
