# Current-centred PDO clinical, AUCrel, and response-extreme associations

## Scope and canonical inputs

This terminal workflow relates baseline malignant PDO composition and activity
to the variables in `PDO_ClinicalMetadata_V3.xlsx`. It uses only the current
centred-refined state vector, current finalized UCell matrix, finalized MP gene
lists, and strict MP order documented in `analysis/ANALYSIS_MAP.md`. It excludes
`SUR843T3_PDO` and all `_Treated_PDO` samples. An `_Untreated_PDO` sample is a
valid baseline sample. The script stops if more than one baseline sample maps
to the same SUR identifier, preventing patient-level pseudoreplication.

The biological and statistical unit is one baseline PDO sample. Cells are
never treated as independent replicates in an association test.

## Clinical variables and normalization

The categorical analyses cover Gender, age (`<=60` versus `>60`, matching the
legacy final clinical workflow), Tumour location, Tumour type, Histology,
AJCC, Clinical response at OG MDT, Mandard tumour regression score, Response
based on Mandard, PDO origin, and pre/post chemotherapy timepoint. Blank and
explicit `N/A` values are excluded per variable. Recorded workbook categories
are otherwise retained. The following display-only harmonizations are made:

- `R`/`Responder` and `NR`/`Non-responder` are normalized consistently;
- PDO origin `T` is displayed as `Primary tumour` and `LN` as `Lymph node`.

For each clinical level, state-composition stacked bars show the mean of the
sample-level proportions, so every PDO receives equal weight. The stacked bars
display the six canonical states plus Hybrid and Unresolved and therefore sum
to 100%. Labels report cell count, sample count, batch count, and the represented
batches. Canonical state percentages use all cells in that sample as their
denominator; only the six canonical state numerators are tested in association
boxplots.

MP activity is the mean finalized UCell score per sample and MP. For two-level
categorical variables the comparison is a two-sided Wilcoxon rank-sum test; for
variables with more than two observed levels it is Kruskal-Wallis. Benjamini-
Hochberg correction is performed across the six states or finalized MPs within
each clinical variable. Raw source data and both raw and adjusted P values are
saved even when groups are too small for useful inference. Following the scRef
final-boxplot convention, large asterisks above features mark raw P below 0.05;
the corresponding BH-adjusted values remain available in the source tables.

## Continuous FLOT AUCrel analysis

The script requires all 12 numeric values in `AUCrel (FLOT)`. For every state
proportion and MP sample mean, it draws the 12 labeled PDO points with an
ordinary least-squares line and 95% confidence band. Fits with raw P below 0.05
are highlighted in red and marked with an asterisk. It saves the slope, 95%
confidence interval, P value, R-squared, and BH-adjusted P value. These plots
describe linear patterns; they do not imply that bounded proportions or the
small cohort support a definitive predictive model.

Higher AUCrel is interpreted as greater FLOT resistance.

## Prespecified response extremes

The response-extreme comparison is fixed before testing:

- sensitive: SUR680, SUR629, SUR1090;
- resistant: SUR727, SUR1181, SUR1363.

These are respectively the three lowest and three highest numeric AUCrel PDOs
in the supplied workbook. Overall state proportions and MP mean UCell scores
are compared with unpaired two-sided Wilcoxon rank-sum tests. Boxplots retain
their conventional median and quartiles, overlay the arithmetic group mean as
a white diamond, and display each PDO as a neutral-coloured point positioned
only over its correct response group. Because each arm contains only three
PDOs, the reported response-extreme effect is the resistant-group arithmetic
mean minus the sensitive-group arithmetic mean. Individual PDO points and
effect direction are more informative than nominal P values; no claim of
population-level clinical inference should be made. The response-extreme
composition stack also includes Hybrid and Unresolved and uses the
equal-PDO-weighted group mean.

State-resolved MP activity is the mean UCell score within each sample-state.
A sample-state is retained when it contains at least 20 cells by default
(`PDO_CLINICAL_MIN_CELLS` can change this auditable threshold). The same
unpaired sample-level test is applied within state and MP, with BH correction
across MPs separately within each state.

## State-resolved selected pathway scoring

For the six response-extreme PDOs, raw RNA counts are summed for every
sample-state pseudobulk. TMM normalization and log2 counts-per-million with a
prior count of 2 are calculated using edgeR. Each gene is z-scored across all
available response-extreme sample-state pseudobulks, and a signature score is
the arithmetic mean of its available member-gene z-scores. Sample-states with
fewer than 20 cells by default are not tested or plotted.

The selected pathways reproduce the matched-response workflow:

- Hallmark E2F targets, G2M checkpoint, apoptosis, p53, DNA repair,
  TNF/NF-kB, EMT, hypoxia, xenobiotic metabolism, unfolded protein response,
  oxidative phosphorylation, and combined interferon-alpha/gamma response;
- `CCSIG`, the 50 most expressed consensus cell-cycle genes in these
  pseudobulks;
- intestinal metaplasia from current MP15 genes;
- ciliated progenitor epithelium from current MP17+ genes.

Within each state and signature, the plotted effect is mean resistant score
minus mean sensitive score. Wilcoxon tests use PDO sample-state scores, and
BH correction is across signatures within each state. The exact gene sets and
pseudobulk logCPM matrix are retained in live storage so figures can be
regenerated without reloading the Seurat object.

## Outputs, cache, and limitations

All source data, statistics, normalized pseudobulk matrices, gene sets,
figures, reports, and logs are written under
`PDOs_outs/Auto_centred_clinical_auc_response_associations/`. The persistent
result cache supports `PDO_REPLOT_ONLY=1`; `PDO_FORCE_REBUILD=1` recomputes the
analysis.

The cohort is small, several categorical levels contain few PDOs, clinical
variables are correlated, and the tests are univariable. Results are
exploratory associations, not adjusted causal or predictive estimates.
