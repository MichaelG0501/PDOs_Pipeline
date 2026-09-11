# Centred high-resolution matched-FLOT PDO metaprograms

## Scope

This terminal matched-treatment workflow deliberately preserves many smaller
metaprograms rather than reducing them to the major canonical PDO states. It is
separate from the canonical centred-refined MP/state route and must not replace
`merged_refined_mp_genes.rds` or the current five-state vector.

The workflow uses the eight explicitly matched samples from SUR1070, SUR1072,
SUR1090 and SUR1181. Each patient contributes one untreated and one FLOT-treated
PDO sample. `SUR843T3_PDO` is excluded by the shared configuration and is not
eligible for this analysis.

## Run order

1. `Auto_pdo_flot_matched_geneNMF.R` runs matched-sample centred GeneNMF.
2. `Auto_pdo_flot_matched_highres_mp_trend_filter.R` extracts diverse MPs,
   scores cells and applies the paired treatment-trend filter.
3. `Auto_pdo_flot_highres_enrichment_annotation.R` annotates retained MPs.
4. `Auto_pdo_flot_highres_cluster_heatmap.R` plots retained MPs directly
   within the current five centred PDO states and visualises their MP-by-MP
   correlation structure.
5. `Auto_pdo_flot_matched_survival_and_state_plots.R` tests the retained MP
   signatures in TCGA EAC survival data.

The two historical per-cell DGE/general response scripts remain legacy. Their
valid current functionality is implemented with paired pseudobulk edgeR and
current states in `centred_refined_pdo_flot_matched_response.R`.

## Centred NMF and high-resolution nMP choice

GeneNMF uses `k=4:9`, `min.exp=0.05` and `center=TRUE`. GeneNMF therefore
centres each gene across cells and sets negative centred values to zero before
factorization, matching the transformation used by the canonical centred PDO
route. The NMF programme object is persistent live input to later steps.

The diverse metaprogram count is fixed by the inherited rule
`round(total NMF programmes / 2)`. This is intentionally not the canonical
rank-selection optimum: it targets approximately two individual NMF
programmes per MP to preserve treatment-responsive heterogeneity. Extraction
uses cosine similarity, specificity weight 5, explained weight 0.5 and minimum
confidence 0.5.

## UCell scoring and paired trend filter

All eight post-QC matched PDO objects are read from `PDOs_list_PDOs.rds`.
Counts are restricted to their common genes, cell identifiers are prefixed by
sample, and MP signatures are scored with UCell (`maxRank=1500`). Mean and
median UCell scores are calculated separately for every MP and sample.

For each patient, an MP receives an `increase` call only when both its sample
mean and sample median are higher after treatment. It receives a `decrease`
call only when both are lower. Mixed mean/median directions do not support
either call. An MP is retained when the same direction is supported by at
least three of four patient pairs. No P-value or effect-size threshold is used
for selection. Paired Wilcoxon values and BH-adjusted values are descriptive
audit fields; four pairs do not provide well-powered population inference.

## Data-derived annotation without manual grouping

Retained MPs are never assigned to manually curated functional groups. Each
MP's short display label is its best 3CA enrichment match after excluding
cell-cycle references. If no eligible overlap exists, the label states that no
non-cell-cycle 3CA match was found. GO, Hallmark, 3CA and developmental
enrichment results are retained in full for interpretation.

The state-resolved heatmaps display individual retained MPs. Rows are split by
the observed increase/decrease selection direction and clustered within those
slices using their plotted score profiles. The columns use the independent,
current centred five-state vector. States do not select, rename or group these
high-resolution MPs.

The MP correlation heatmap is computed independently within each of the eight
matched patient-treatment samples using Spearman correlation between per-cell
UCell scores. Correlations are clipped only for the Fisher transformation,
combined across samples by their mean Fisher-Z value, and transformed back to
Spearman rho. A one-sample t-test of the finite sample-level Fisher-Z values
against zero is retained as a descriptive p-value matrix. MPs are clustered
across the complete correlation matrix using average linkage on `1-rho`;
treatment direction is an annotation only and does not split or constrain the
clustering. This sample-balanced construction prevents samples with more cells
from dominating the displayed MP relationships.

## TCGA survival analysis

Selected MP gene sets are intersected with reconstructed TCGA ESCA TPM and
must retain at least five genes. GSVA scores are calculated for primary EAC
cases. Univariate Cox associations are reported for continuous scores, median
splits and upper-versus-lower quartile splits, following the sensitivity
layout of `Auto_08_tcga_mp_survival_volcano_centred.R`. BH adjustment is
performed within split method. These are outcome associations, not evidence
that the MPs predict FLOT response or cause survival differences.

## Storage and cache policy

Programme objects, scores, cell metadata, retained genes, enrichment, source
tables and survival results are persistent under
`PDOs_outs/Auto_pdo_flot_centred_highres_metaprogram_trends/`. They are
sufficient to regenerate downstream figures without repeating NMF. Analytical
caches are reused by default and rebuilt with `--force` or
`PDO_FORCE_REBUILD=1` where documented. All analytical execution uses PBS.
