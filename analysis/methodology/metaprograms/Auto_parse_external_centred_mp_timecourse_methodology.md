####################
# PDO/scRef centred-refined MP activity across Parse timepoints

## Status and scope

This is an active terminal cross-dataset comparison. It projects the current
centred-refined MP gene sets from the PDO and scRef pipelines into the Parse
NACT1090 time course. It does not derive, rename, merge, or match MPs.

## Current MP definitions

PDO signatures and plotting order come from:

- `PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds`
- `PDOs_outs/centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv`
- `analysis/metaprograms/centred/Auto_05_centred_refined_mp_ordered_heatmaps.R`

scRef signatures and plotting order come from:

- `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`
- `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/tables/centred_refined_mp_state_grouping.csv`
- `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/analysis/metaprograms/centred/05_centred_refined_mp_ordered_heatmaps.R`

MP numbers have no cross-pipeline meaning. Every score column is therefore
source-prefixed (`PDO__<MP>` or `scRef__<MP>`), and every plot label includes
the source. The two scRef rows whose current state is explicitly `Excluded`
(`MP11c`, `MP18a`) are documented but omitted from active scoring.

## Parse inputs and scoring

The six timepoints are shown in the fixed order `T0`, `T1`, `T2`, `T4`, `R4`,
`eR4`. The input for each is the live final Seurat object:

`/rds/general/project/spatialtranscriptomics/live/Parse_Pipeline/parse_outs/by_samples/<sample>/Auto_<sample>_final.rds`

RNA counts are scored independently per timepoint using
`UCell::ScoreSignatures_UCell`, `maxRank = 1500`, and six cores. A signature
must have at least five detectable genes in every sample. Per-sample caches,
the combined per-cell score matrix, exact signature lists, and aligned cell
metadata are persisted under the live PDO project so figures remain
reproducible without ephemeral storage.

## Biological ordering and displays

Rows are not clustered. They are ordered by:

1. Cell cycle
2. Classic proliferation
3. Basal to intestinal metaplasia
4. SMG to intestinal metaplasia
5. Stress adaptive
6. Cancer-cell immune mimicry
7. PDO medium induced

Within a biological group, PDO MPs precede scRef MPs while retaining the
current order defined by each source workflow. Heatmaps contain both an MP
source annotation and a biological-group annotation.

The report contains:

- a within-MP temporal z-score heatmap, appropriate for comparing time-course
  shapes between signatures;
- a raw mean-UCell heatmap, which retains absolute UCell score scale;
- one mean/IQR trend page per biological group;
- one median/IQR and 5th-95th percentile distribution page per group.

UCell scores are enrichment/activity scores, not direct expression values.
Raw score magnitude can depend on signature composition, so cross-signature
interpretation should emphasize temporal patterns and not treat absolute
differences as equivalent gene-expression fold changes.

## Statistical summaries

For every MP, the workflow reports a six-timepoint Kruskal-Wallis test and a
two-sided Wilcoxon comparison of pooled `T2+T4` versus pooled `T0+eR4`.
The larger pooled group is deterministically downsampled to the size of the
smaller group. P values are Benjamini-Hochberg adjusted separately for each
test family. These cell-level tests describe this dataset and do not substitute
for biological-replicate inference.

## Cache controls and outputs

`PDO_FORCE_REBUILD=1` recomputes all UCell caches.
`PDO_REPLOT_ONLY=1` reuses the persistent combined score and metadata caches.

All outputs are under:

`PDOs_outs/Auto_parse_external_centred_mp_timecourse/`

The combined signatures, per-cell UCell matrix, metadata, and plot matrices
are reusable downstream inputs. PDFs and PNGs are terminal presentation
figures.
####################
