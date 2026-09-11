# PDO–scRef centred-state concordance methodology

## Scientific question

This analysis asks whether transcriptional states defined independently in PDO
and scRef/scATLAS data reproduce one another across datasets. It does not
project old scRef MP numbers onto PDO cells and does not use the removed
`PDOs_final.rds` route.

The current step-06 named state vectors are treated as fixed. The current
step-05 MP-to-state tables and merged refined MP gene lists define independent
PDO and scRef state signatures. MP numbers are dataset-local identifiers and
are never matched between datasets.

Four states have an a priori one-to-one correspondence:

- Classic proliferation
- Basal to intestinal metaplasia
- SMG to intestinal metaplasia
- Stress adaptive

`PDO medium induced` and `Cancer-cell immune mimicry` are dataset-specific and
are retained as negative-control alternatives rather than forced into a fifth
match.

## Evidence layers

### Current MP signature cross-projection

Genes from all current centred refined MPs assigned to a state are combined
within each dataset. Both the PDO-derived and scRef-derived state signatures
are scored in both datasets with UCell. Mean scores are calculated within the
independently assigned receiving states and standardized within each signature
across receiving states. The two directional heatmaps therefore test whether a
state signature derived in one dataset peaks in the independently named state
of the other dataset.

The same current state signatures are also compared directly at gene level.
Every PDO/scRef state pair receives a Jaccard overlap and one-sided Fisher exact
test over the union of genes in the current state signatures. This provides a
score-free check unaffected by UCell scaling or MP identifiers.

### Current scRef state signatures projected onto PDO cells

Every current scRef MP assigned to one of the five scRef step-06 states is
scored separately in PDO cells. The resulting MP scores are normalized using
the same sample-centering, study-scaling, maximum-within-state, threshold
(`0.5`), and hybrid-gap (`0.3`) rules as the scRef step-06 state definition.
The projected scRef-signature state is compared with the independently assigned
PDO step-06 state in a row-normalized stacked barplot. Projected `Unresolved`
and `Hybrid` cells remain visible, preventing uncertain assignments from being
removed before the matched percentage is calculated.

### Sample-aware pseudobulk state effects

Raw counts are summed within each biological sample and state, retaining groups
with at least 20 cells. State-vs-rest effects are fitted independently within
each dataset with edgeR normalization and limma-voom models containing sample
fixed effects. This makes samples—not individual cells—the inferential unit and
controls sample-level baseline differences.

The cross-dataset matrix contains Spearman correlations between independently
estimated state-vs-rest log-fold-change vectors over a fixed union of
informative genes. This is more resistant to global platform and tissue
baseline shifts than directly correlating raw state-average expression.

### Independent marker overlap

For each state and dataset, the top 200 positive genes ranked by moderated
t-statistic are selected independently. Every PDO/scRef state pair is assessed
with Jaccard overlap and a one-sided Fisher exact test over the genes tested in
both datasets. The full matrix is retained; expected matches are not selected
after seeing the results.

### Specificity test

For the four prespecified shared states, the mean expected diagonal score is
compared with all alternative cross-state pairings. Exact enumeration of all
24 scRef state-label permutations supplies a mapping-level permutation
p-value. Reciprocal best matches, the diagonal margin, and the unmatched
dataset-specific states are reported explicitly.

For each shared state and evidence layer, the expected-match value minus the
best alternative-state value is also displayed. Positive values support a
state-specific match; negative values show that another state is at least as
similar. This prevents a strong overall diagonal statistic from concealing an
unsupported individual state.

The state-vs-rest transcriptional effects are additionally displayed as five
pages, one per PDO state. Each page compares that PDO state with all five scRef
states in separate panels using common axes and genes. This exposes competing
off-diagonal correlations rather than showing only the prespecified matched
pairs.

## Interpretation

No single matrix proves identity. Evidence for state conservation is strongest
when the same prespecified correspondence is supported by:

1. reciprocal current-MP signature projection;
2. direct overlap of current MP-state genes;
3. state-vs-rest transcriptomic effect correlation;
4. independently ranked marker overlap; and
5. a positive diagonal-specificity margin under exact label permutation.

The design follows the cross-dataset prediction principle of MetaNeighbor,
which defines replicability as the ability of transcriptional identity learned
in one dataset to recover the corresponding population in another:
https://www.nature.com/articles/s41467-018-03282-0

Sample-level pseudobulk inference is used because benchmarking shows that
aggregation within biological replicates controls false discoveries better
than treating cells as independent:
https://www.nature.com/articles/s41467-021-25960-2

## Outputs and replotting

The presentation PDF, all matrices, full DGE results, overlap statistics,
cell-count audit, and run log are saved under
`PDOs_outs/Auto_PDO_scRef_centred_state_concordance/`. The live RDS cache
contains the pseudobulk inputs and derived matrices required to regenerate the
figures. Set `PDO_REPLOT_ONLY=1` to regenerate tables and figures from this
cache, or `PDO_FORCE_REBUILD=1` to ignore it.
