# Centred-refined matched-FLOT PDO response

This workflow replaces the legacy nMP13 matched-FLOT visualization with the finalized centred-refined PDO MP model and its noreg state assignments. It retains the four explicitly matched untreated/FLOT PDO pairs (SUR1070, SUR1090, SUR1072 and SUR1181); no unmatched PDO is introduced into a paired comparison.

State composition is calculated within each sample. The publication panel reports paired percentage-point changes, with a paired Wilcoxon effect summary across the four PDO pairs. Because four pairs cannot yield a two-sided exact Wilcoxon P value below 0.05, the figure emphasizes patient-level effects and medians and does not use significance stars.

MP effects are within-patient changes in mean UCell score for the strict finalized MP order. The workflow does not read the legacy nMP13 GeneNMF object, `UCell_scores_filtered.rds`, `Auto_PDO_mp_adj_noreg.rds`, or `Auto_PDO_final_states.rds`.

For pathway response, raw RNA counts are aggregated by sample and finalized state. A pseudobulk is retained when it contains at least 20 cells by default. TMM-normalized logCPM values are gene-wise z-scored, and Hallmark scores are means across available member genes. Treated-minus-untreated pathway deltas are paired within patient and state. The exact Hallmark gene sets used and the pseudobulk matrix are saved persistently under `PDOs_outs/centred_refined_flot_matched_response/intermediate/` so all panels can be replotted without reloading the full PDO object.

State-resolved differential expression uses paired edgeR quasi-likelihood models (`~ Patient + Treatment`) and requires at least three complete patient pairs. All result tables are retained regardless of significance. These exploratory molecular response results should not be described as population-level clinical inference because the biological unit is four PDO pairs.

Figures are written as editable vector PDFs at Nature double-column width (183 mm), using Arial, explicit 6.5–7 pt typography, color-blind-aware state colours and a 600-dpi composite PNG. Exact panel source-data CSV files are written under the live output directory.
