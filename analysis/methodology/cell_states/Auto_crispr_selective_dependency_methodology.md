# Oesophageal CRISPR Selective-Dependency Analysis for PDO MPs and States

## Scope and interpretation

This terminal workflow evaluates whether genes defining the canonical centred
PDO metaprograms (MPs) and states overlap the significant oesophageal
differential-dependency genes reported by Herranz-Ors et al. (2026;
https://doi.org/10.1038/s41586-026-10830-y). It uses
only the `diff_dep_analysis_oesophageal` sheet of
`CRISPR_screen_results_PDO.xlsx`; colorectal and pooled gastrointestinal sheets
are not analysed.

The workbook is a filtered context-specific dependency resource, not a complete
essentiality screen. The source study removed pan-cancer and organoid core
fitness genes, reference essential/non-essential controls, non-expressed genes,
and genes with nearly uniform depletion status before differential-dependency
testing. Consequently:

- overlap supports a selective OAC knockout vulnerability among genes that
  define an MP or state;
- strong, recurrent overlap suggests that the transcriptional identity may be
  selectively vulnerable rather than universally essential;
- absence means no selective-dependency evidence in this filtered workbook and
  is interpreted operationally as limited selective lethal targetability;
- absence does not prove biological non-essentiality, because a gene may have
  been excluded as uniformly essential, uniformly non-depleted, non-expressed,
  or otherwise outside the differential-dependency analysis.

This workflow does not infer that cells currently occupying a state will be
depleted in a pooled CRISPR experiment. It prioritizes state/MP programs whose
definition genes contain knockout dependencies measured across independent
oesophageal organoids.

The same analysis is performed independently for the current scATLAS MPs and
states. scATLAS signatures are assessed against a scATLAS-derived expression
universe rather than the PDO universe, so their expression-matched null results
are dataset-specific.

## Canonical PDO inputs

The analysis uses the current centred route only:

- `PDOs_outs/centred_mp_refinement/merged_refined_mp_genes.rds`;
- `PDOs_outs/centred_mp_refinement/merged_refined_mp_gene_weights.rds`;
- `PDOs_outs/centred_mp_refinement/centred_refined_noreg_states.rds`;
- `PDOs_outs/PDOs_merged.rds`;
- MP/state grouping, order, descriptions, and exclusions from
  `analysis/shared/Auto_pdo_analysis_config.R`.

`SUR843T3_PDO`, `Hybrid`, and `Unresolved` cells are excluded. Cell-cycle MPs
are reported individually for QC and biological context but do not contribute
to combined state gene lists.

## Oesophageal CRISPR reference

All 4,082 rows in the oesophageal workbook sheet have FDR-adjusted P values
below 0.05. Gene identifiers are matched first by case-insensitive HGNC symbol
and then, when applicable, by version-stripped Ensembl identifier.

For every matched gene, the analysis retains all supplied fields and derives:

- depletion prevalence = `N_depleted / (N_depleted + N_not_depleted)`;
- conditional depletion strength = `max(0, -MEAN_LFC_depleted_group)`;
- prevalence-weighted strength = depletion prevalence multiplied by conditional
  depletion strength;
- overall depletion strength = `max(0, -MEAN_LFC)`.

The prevalence-weighted strength is the primary gene-level depletion burden.
It is high only when knockout is both frequent across the 59 oesophageal
organoids and strong within the depleted group. Overall mean LFC burden is
retained as a sensitivity measure.

## Canonical scATLAS inputs and signatures

The independent scATLAS analysis uses `EAC_Ref_epi.rds`, the current centred
refined MP genes and weights, the current centred refined noreg state vector,
the MP/state groupings in `analysis/shared/scRef_config.R`, and the ranked state
markers produced by `analysis/cell_states/final_state_marker_discovery.R`.
Only cells assigned to the five canonical biological states contribute to the
scATLAS expression universe; Hybrid and Unresolved cells are excluded.

All 17 current scATLAS MPs are scored individually. For each of the five states,
the genes of all component MPs are also combined and deduplicated. A second
state definition uses the top 150 genes in the current ranked marker table,
ordered exactly as in the marker-discovery workflow. Positive median
within-sample log2 fold change provides the definition weight for these marker
signatures.

The scATLAS expression universe includes genes detected in at least 10
canonical-state cells. Expression bins, matched permutations, overlap metrics,
depletion-burden metrics, empirical P values, FDR correction, and evidence
classes are then calculated exactly as for PDO, but independently within this
scATLAS universe.

## Signature definitions

### Individual MPs

Each of the 15 canonical centred MPs uses its complete step-04 gene list. NMF
weights from `merged_refined_mp_gene_weights.rds` are normalized to sum to one
within each MP and provide the primary definition weights. Unweighted overlap
fractions are also reported.

### Combined state-component MP genes

For each of the six biological states, the complete genes from all non-cell-
cycle MPs assigned to that state are combined and deduplicated. A gene occurring
in more than one component MP is counted once. All combined genes receive equal
weight so that the result represents the state gene-list composition rather
than the number or scale of its component MPs.

### Fresh state edgeR signatures

The current state vector is used to construct a fresh state-vs-rest signature
for each of the six states. Raw RNA counts are aggregated into target and rest
pseudobulks separately within each PDO sample. A pseudobulk must contain at
least 20 cells, and a sample is retained only when both target and rest pass
that threshold. At least three paired samples are required.

For each state, edgeR uses:

1. `filterByExpr()` with the sample-blocked design;
2. TMM normalization;
3. robust dispersion estimation;
4. robust quasi-likelihood fitting;
5. the coefficient for target state versus pooled other canonical states in
   the design `~ sample + group`.

This mirrors the sample-blocked pseudobulk strategy in
`Auto_drug_reversal_inputs.R` but uses the current six centred states. The
state-defining DGE list is the top 150 positive-logFC genes, ordered by FDR,
then decreasing logFC and target-cell detection. This deliberately reproduces
the drug-reversal signature convention. Positive logFC values, normalized to
sum to one within a state, provide definition weights. Full, untruncated edgeR
results remain persistent live outputs.

## Classic-proliferation cross-dataset three-way overlap

An additional four-panel analysis compares current scATLAS, oesophageal CRISPR,
and current PDO Classic-proliferation gene sets. The scATLAS inputs are the
centred refined `MP2+` gene list and the top 150 Classic-proliferation genes from
`ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv`.
The latter is produced by
`scRef_Pipeline/analysis/cell_states/final_state_marker_discovery.R` and retains
that workflow's ordering by decreasing ranking score, reproducibility, median
within-sample log2 fold change, and specificity gap. The PDO inputs are the
centred refined `MP19+` gene list and the top 150 Classic-proliferation genes
from the sample-blocked edgeR workflow described above.

Classic proliferation is defined by a single MP in each atlas (`MP2+` in
scATLAS and `MP19+` in PDO), so a separate combined-state MP set would be an
exact duplicate and is not plotted. The four panels are the Cartesian product
of the two scATLAS definitions and two PDO definitions, with the 4,082-gene
oesophageal CRISPR list fixed as the third set. Gene identifiers are
case-insensitive; genes matching the CRISPR reference by symbol or
version-stripped Ensembl identifier are canonicalized to the CRISPR HGNC
symbol.

Each panel reports all seven mutually exclusive regions of a three-set overlap,
including the triple intersection. Circles are deliberately schematic and not
area-proportional. The exact region counts, per-gene memberships, selected gene
sets, input marker ranks, and a dedicated four-panel table of the triple-overlap
genes are saved as source tables/intermediate data, so the figure is auditable
without extracting values from the PDF.

## Overlap and depletion-strength scoring

For every MP or state signature, the workflow reports:

- number and fraction of definition genes in the filtered CRISPR list;
- definition-weighted overlap fraction;
- median depletion prevalence among overlapping genes;
- median depleted-group strength among overlapping genes;
- prevalence-weighted depletion burden across the complete signature;
- overall mean-LFC burden;
- Fisher overlap odds ratio as a simple expression-universe comparison.

Genes absent from the filtered CRISPR list contribute zero to overlap and
depletion burden. They are labelled as lacking evidence in this filtered
resource, not as proven non-essential genes.

## Expression-matched null model

MP and DGE genes are biased toward detectable, often abundant transcripts.
Because screen inclusion is also expression-dependent, raw overlap alone can
be misleading. The analysis therefore builds a PDO expression universe from
genes detected in at least 10 canonical-state cells, calculates overall log
CPM, and divides genes into ten expression bins.

For each signature, 10,000 random gene sets are sampled while preserving its
gene count, expression-bin composition, and definition weights. The observed
weighted overlap and prevalence-weighted burden are compared with these null
distributions. The workflow reports empirical upper-tail P values, z-scores,
Benjamini-Hochberg FDR values within each signature type, and a dependency
index from 0 to 100 equal to the empirical percentile of the primary burden.

The presentation evidence classes use empirical FDR 0.10:

- **Strong**: both overlap and burden are enriched;
- **Depletion-strength enriched**: burden alone is enriched;
- **Overlap enriched**: overlap alone is enriched;
- **Nominal overlap and/or burden enrichment**: the corresponding empirical
  P value is below 0.05 but does not pass FDR correction;
- **Above-background trend**: burden z-score is non-negative but neither the
  nominal nor FDR threshold is met;
- **Below expression-matched background**: burden z-score is negative and no
  enrichment threshold is met.

Point size in the MP/state summary figures represents raw gene-overlap
fraction. Horizontal position represents the expression-matched percentile of
prevalence-weighted depletion burden. Colour represents the evidence class.
These are deliberately separate encodings: a signature can have substantial
raw overlap without reaching multiple-testing-adjusted enrichment.

The continuous components and exact gene-level tables are authoritative; the
class is a presentation summary and must not be read as a clinical prediction.

## Storage, cache, and execution

Full DGE, gene-level matches, score components, null distributions, source
tables, figures, and an Excel review workbook are stored under live
`PDOs_outs/Auto_crispr_selective_dependency/`. They are sufficient to recreate
the presentation figures without the ephemeral filesystem.

Set `PDO_FORCE_REBUILD=1` to rebuild edgeR and scoring caches. Set
`PDO_REPLOT_ONLY=1` to require and reuse the completed scoring cache. Run with
`qsub Auto_run_crispr_selective_dependency.sh` in the `dmtcp` environment.
