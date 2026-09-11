####################
# PDO CNA Subclone Methodology

This workflow assigns CNA subclones within each PDO sample from the PDO-only InferCNA matrix generated with the Carroll 2023 non-malignant reference.

## Inputs

- InferCNA target matrix: `PDOs_outs/cnv/Auto_PDO_infercna_target_outs_Carroll_2023.rds`
- InferCNA target metadata: `PDOs_outs/cnv/Auto_PDO_infercna_target_meta_Carroll_2023.rds`
- Gene order: `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`
- Merged PDO object: `PDOs_outs/PDOs_merged.rds`
- Per-sample PDO objects: `PDOs_outs/by_samples/<sample>/<sample>.rds`
- Final PDO state vector: `PDOs_outs/Auto_PDO_final_states.rds`
- PDO metaprogram scores: `PDOs_outs/Auto_PDO_mp_adj_noreg.rds`

## Cell Matching

InferCNA target columns retain the InferCNA cell identifiers and an `original_cell` barcode in metadata. The downstream PDO state and metaprogram matrices use merged PDO cell names. The script reconstructs those merged names as:

```r
paste(sample, original_cell, sep = "_")
```

when the original barcode is not already present in the merged object. Only cells present in the CNA matrix, merged PDO metadata, final state vector, and metaprogram score matrix are analysed.

## CNA Matrix Preparation

For each sample:

1. Genes are intersected with the hg38 gene-order table.
2. Genes are ordered by chromosome and genomic start.
3. Non-finite rows are removed.
4. Low-signal genes are removed by retaining the upper two thirds of genes ranked by mean absolute InferCNA signal.
5. Genes are binned within chromosomes in consecutive 100-gene bins for clustering and plotting.

## Candidate Subclone Generation

The filtered CNA matrix is transposed to cells by genes and clustered with Louvain community detection on a 15-nearest-neighbour graph (`k = 15`). No PCA/Ward candidate cuts and no silhouette filter are used for subclone definition.

The sample-level analysis minimum is:

- minimum sample cells for analysis: 40

After Louvain clustering, undersized provisional clusters are dropped from the subclone analysis rather than merged into a larger clone. A provisional cluster must have at least:

- 20 cells
- 5% of the sample cells

The larger of these two thresholds is used for each sample. This is intentionally relaxed so small subclones can be retained when there are enough cells to see a coherent CNA pattern. Dropping undersized clusters prevents very small Louvain communities from becoming false-positive subclones while also preventing their CNA signal from being absorbed into true larger subclones.

The remaining Louvain clusters are treated as provisional CNA clusters and are then merged by chromosome-arm CNA patterns.

## Distinctness Checks

For every provisional cluster, the script calculates mean CNA values across genes in each chromosome arm. A chromosome arm is called:

- amplified if the cluster mean is > 0.10
- deleted if the cluster mean is < -0.10
- neutral otherwise

This less brittle arm-call threshold is used for the PDO and SCREF reruns because several visually clear subclones had arm means below ±0.15 and were therefore collapsed to neutral/neutral patterns.

The script then iteratively compares all cluster pairs. A pair is merged only when there is no robust chromosome-arm evidence that they are distinct. Distinct evidence is defined as either:

- at least one chromosome arm has different amplified/deleted/neutral calls and an absolute mean CNA difference of at least 0.08
- at least two chromosome arms have absolute mean CNA differences of at least 0.08

Pairs are merged when their maximum chromosome-arm mean difference is below 0.06. Pairs are also merged when they show the same CNA shape with only a strength difference: centred chromosome-arm profile correlation at least 0.80 and mean absolute arm difference no more than 0.015. This prevents global CNA signal-strength variation, such as in `SUR680T3_PDO`, from being miscalled as separate subclones while retaining pattern-shifted samples such as `Carroll_2023_EAC-ACMO_ICI-4W_tumour_frozen`.

After each merge, arm means and calls are recomputed, and the pairwise comparison repeats until no further mergeable pair remains.

## Selecting The Final Number Of Subclones

The remaining merged Louvain clusters are the final CNA subclones. There is no silhouette-based selection step.

This follows the published Nature-paper subclone approach provided by the user:

- retain the top 67% CNA-signal genes
- generate Louvain clusters with `k = 15`
- drop undersized Louvain clusters instead of merging them into larger clones
- call chromosome-arm amplification/deletion/neutrality at ±0.10
- iteratively merge clusters without robust chromosome-arm CNA evidence
- test the remaining subclones for MP-expression differences

The per-sample selected solution and all candidate diagnostics are written to:

- `PDOs_outs/Auto_PDO_cnv_subclone_mp/Auto_PDO_cnv_subclone_summary.csv`
- `PDOs_outs/Auto_PDO_cnv_subclone_mp/Auto_PDO_cnv_subclone_cluster_diagnostics.csv`

## QC Metrics

RNA QC metrics are read from the merged PDO object when present. `percent.mt` was not retained in `PDOs_merged.rds`, so the script restores it from each per-sample PDO object in `PDOs_outs/by_samples/<sample>/<sample>.rds`, where the original QC pipeline stored `percent.mt`.

## Plot Colour Consistency

Subclone colours are fixed by label and reused in every plot:

- Subclone A: red
- Subclone B: blue
- Subclone C: green
- Subclone D: purple
- Subclone E: orange
- Subclone F: brown

The same mapping is used in CNA heatmap annotations, metaprogram boxplots, metaprogram mean heatmaps, and QC/CNA boxplots.
####################
