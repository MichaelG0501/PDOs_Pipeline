# Centred PDO Sample-Abundance Methodology

`sample_abundance_pdo.R` is a terminal plotting script. It aligns the named
canonical centred state vector and the current centred-refined UCell matrix to
`PDOs_merged.rds` by exact cell barcode, then summarizes cells by PDO sample.

Samples are displayed with clinical annotations read from the live clinical
workbook. Batch labels are derived from sample names: treated, untreated,
new-four-sample, or historical PDO. Age is split at 60 years (`<=60`, `>60`);
missing or unrecognized clinical values are displayed as `Unknown`. These are
display categories, not fitted statistical thresholds.

State order, MP descriptions, and colors follow the canonical centred ordering.
The state panel shows per-sample proportions including `Unresolved` and
`Hybrid`; the MP panel summarizes current merged-refined UCell activity. Cell
counts are descriptive and are not treated as independent patient-level
replicates. The PDF is terminal and no output from this script is a downstream
state or score input.
