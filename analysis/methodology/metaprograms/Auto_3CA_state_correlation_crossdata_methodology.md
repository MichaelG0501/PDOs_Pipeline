# State-resolved cross-dataset 3CA correlation methodology

## Aim

This terminal analysis compares the five centred, refined, non-regressed PDO
states with the five centred, refined, non-regressed scRef/scATLAS states using
the same per-cell 3CA UCell score features in both datasets.

## Inputs and inclusion rules

The analysis reads the live `UCell_3CA_MPs.rds` score matrix and named
`centred_refined_noreg_states.rds` vector from each pipeline. Cell barcodes are
matched within each dataset. Only the five target states explicitly defined in
the corresponding state-definition script are included; `Unresolved` and
`Hybrid` cells are excluded. No UCell rescoring or state reassignment is
performed.

## State-level comparison

For every target state, the arithmetic mean UCell score is calculated for each
3CA metaprogram shared by the PDO and scRef score matrices. Every PDO-state and
scRef-state pair is compared by Spearman correlation across these shared 3CA
metaprogram mean-score vectors.

The first PDF page is a 5 by 5 matrix of the resulting Spearman rho values.
Pages 2 through 6 correspond to successive PDO states. Each page contains five
scatter panels, one against each scRef state in the state-definition order.
Each point is one shared 3CA metaprogram. The axes are synchronized across all
25 panels, dotted lines show the original cross-data score threshold of 0.1,
and metaprograms reaching 0.1 in either member of a state pair are labelled.

## Outputs

The single multi-page presentation PDF is written under
`PDOs_outs/Auto_3CA_state_correlation_crossdata/figures/`. The state mean
scores, 5 by 5 matrix, pairwise statistics, and complete scatter plotting data
are written under `tables/`, allowing the figure to be audited or regenerated
without recomputing state means. A package/session run summary is written under
`logs/`. All outputs are terminal and are not current inputs to another script.
