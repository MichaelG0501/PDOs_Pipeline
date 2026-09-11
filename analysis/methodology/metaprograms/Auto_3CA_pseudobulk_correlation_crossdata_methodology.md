# 3CA Pseudobulk Cross-Dataset Correlation Methodology

## Comparison design

The workflow compares external pan-cancer 3CA signatures across three data
types: scATLAS EAC epithelial cells, OAC PDO single-cell data aggregated by PDO
sample, and OSCC/ESCC tumour-organoid bulk TPM from GSE269447. It does not use
PDO GeneNMF labels and is therefore independent of the centred-vs-uncentred PDO
MP choice.

## Gene sets and score construction

3CA signatures are read from the live `New_NMFs.csv`. Empty and missing entries
are removed. A signature is scored in a dataset only when at least five genes
are present, preventing near-empty overlaps from being interpreted as the
original program. The overlap table records the available gene count for every
signature and dataset.

Single-cell datasets are scored with rank-based UCell and aggregated to a mean
score per biological sample. For bulk RNA-seq, the script constructs the
comparable signature summaries implemented in its bulk-scoring function after
gene-symbol harmonization. Dataset comparison uses the mean across biological
samples, not a pooled cell-weighted mean.

## Correlation and display thresholds

Spearman correlation is calculated across 3CA MPs for each reference-target
dataset pair, with `exact=FALSE` because tied rank scores are expected. A point
is labelled/highlighted when either axis has mean score at least 0.1. This 0.1
cutoff is a readability rule for label density; it is not a statistical or
biological activity threshold. All MPs, including lower-score points, remain in
the correlation and exported tables.

Pan-cancer breadth categories come from the supplied coverage summary:
`General` for more than 12 cancer types, `Shared` for 6–12, and `Specific` for
at most 5. These categories annotate points and do not alter scores or fits.

## Interpretation and storage

Correlation measures concordance of relative 3CA program activity across MPs;
it does not establish cell-state equivalence or treatment effects. Histology,
platform, and bulk-versus-single-cell composition can all contribute to
differences. All score, overlap, and summary tables plus figures are saved in
live storage so the comparison is auditable without relying on ephemeral data.
