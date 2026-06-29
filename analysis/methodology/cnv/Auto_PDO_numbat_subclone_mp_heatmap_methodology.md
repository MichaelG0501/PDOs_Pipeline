# Auto_PDO_numbat_subclone_mp_heatmap methodology

Status: active optional CNV validation, conservative-focused for current PDO subclone analysis.

This workflow combines Numbat clone calls with PDO metaprogram and final-state activity. It reads the Numbat manifest and per-sample clone/cell-map/joint posterior files, projects Numbat posterior CNV calls onto genome-wide bins for per-sample heatmaps, and summarizes MP/state heterogeneity across Numbat subclones.

Current significance calling for MP expression differences uses pairwise subclone comparisons. For each sample, each metaprogram, and each displayed subclone pair, the effect score is:

`median(subclone_1 MP score) - median(subclone_2 MP score)`, divided by the sample-level MAD of that MP score.

The cohort-level diagnostic distribution is built separately for each metaprogram from absolute pairwise scores. The current threshold is `max(1, median(abs(score)) + 3 * MAD(abs(score)))`, with a q95 fallback if the robust estimate is unavailable. A pairwise subclone MP difference is called significant only when it passes both the MP-specific score threshold and Wilcoxon BH FDR < 0.05.

Primary outputs are written under `PDOs_outs/Auto_PDO_numbat_subclone_mp_conservative/` when `PDO_NUMBAT_CLONE_MODE=conservative`, including the sample-page PDF, the two-page cohort summary PDF, the MP-specific diagnostic PDF, pairwise test tables, and threshold tables.
