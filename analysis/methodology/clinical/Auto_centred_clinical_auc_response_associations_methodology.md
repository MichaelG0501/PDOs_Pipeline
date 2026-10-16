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

## Standard pseudobulk differential expression and pathway enrichment extension

The clinical entry point now calls
`analysis/clinical/Auto_centred_clinical_pseudobulk_de.R`. Its dedicated PBS
wrapper, `analysis/clinical/Auto_centred_clinical_pseudobulk_de.sh`, sets
`PDO_CLINICAL_DE_ONLY=1` and executes only this extension. A normal invocation
of the clinical R script executes the previous analyses and then this extension.
This avoids rerunning clinical plots and UMAP when only DE is requested.

The overall comparison aggregates uncorrected raw RNA counts across **all**
baseline cells within each of the six prespecified PDOs, including Hybrid and
Unresolved. The biological replicate remains the PDO/patient, with exactly
three PDOs per group; cells never increase inferential sample size. The script
checks that the fixed groups still match the lowest and highest three AUCrel
values among the 12 independent baseline PDOs. It rejects repeated baseline
PDOs per patient, missing batch metadata, and noninteger or negative
counts. Only Seurat RNA `layer="counts"` is used; no integrated expression,
normalized counts, alternative method, or failed-model fallback is allowed.

The recorded Batch variable has three levels in these six PDOs: `PDO`
(SUR629 and SUR680 sensitive; SUR727 resistant), `Untreated_PDO` (SUR1090
sensitive; SUR1181 resistant), and `new4samples` (SUR1363 resistant only).
The primary fixed-effects design is `~ batch + sensitivity_group`, with
Sensitive as the response reference. This is full rank with six observations,
four coefficients, and two residual degrees of freedom. SUR1363's singleton
batch does not supply a within-batch response comparison. Adjustment applies
only to recorded acquisition groups; it cannot eliminate unrecorded technical
variation or separate effects that are completely confounded. Batch adjustment
is performed in the model; expression matrices are not corrected beforehand.

For each estimable comparison, `edgeR::filterByExpr` uses the actual design,
followed by TMM normalization, `estimateDisp(..., robust=TRUE)`,
`glmQLFit(..., robust=TRUE)`, and `glmQLFTest`. The contrast is Resistant minus
Sensitive: positive log2 fold change means higher expression in resistant
PDOs. Benjamini-Hochberg FDR <0.05 across tested genes within the comparison
is the significance criterion. No fold-change cutoff is imposed. Full
statistics are saved. Only the Overall scope is tested.

The DE extension performs only the overall all-cell comparison. No per-state
pseudobulks or per-state differential expression are calculated. Overall DE
may reflect both state-composition differences and changes in expression
within cells. Canonical state assignments are not required by this extension.

All 50 human MSigDB Hallmark sets from installed `msigdbr` 7.5.1 are tested
using the DGEList method of `limma::camera`, with the same design and contrast,
`inter.gene.cor=0.01`, and at least 10 retained genes per pathway. This
competitive gene-set test uses the expressed gene universe and accounts for
inter-gene correlation. It does not select genes by DE significance before
enrichment. Up indicates a coordinated shift towards resistant PDOs, Down
towards sensitive PDOs. BH FDR is calculated across all 50 Hallmark pathways. Expression enrichment supports an association with pathway activity;
it does not establish biochemical activation or causation. No alternative
pathway method is attempted when dependencies or inputs fail.

All new results live under
`PDOs_outs/Auto_centred_clinical_auc_response_associations/Auto_sensitive_resistant_DE/`.
The **single** figure artifact is
`figures/Auto_sensitive_resistant_pseudobulk_DE_pathways.pdf`. It includes
exactly two figure-only pages: (1) sample MDS and cell counts side by side,
(2) overall batch-adjusted DE volcano and directional Hallmark lollipop plots
side by side. The pathway plot showing only pathways with BH FDR <0.05. Point size represents
the number of tested member genes. All 50 pathways contribute to multiplicity
correction and remain in the source table. If no pathway is significant, the
pathway figure states that result rather than displaying nominal hits.
Legends use separate rows and outer margins to fit within each panel.
There are no text-only pages, state-specific figures, gene heatmaps, repeated
pathway panels, or signature-expression heatmaps. Methods, sample metadata,
model diagnostics, and summary counts remain in the report/source files.
The MDS uses uncorrected TMM logCPM and displays batch by point shape;
inferential DE and CAMERA results remain batch-adjusted.

Raw overall pseudobulks, normalized DGE objects, logCPM, model, gene sets, metadata,
input fingerprints, and summary statistics are cached in live storage. CSV
source data, gene-filter decisions, design matrices, model audits, and logs
also remain in live storage. Replot-only runs require a cache matching input
paths, file sizes and modification times, groups, cell threshold, method
version, and relevant package versions. `PDO_FORCE_REBUILD=1` rebuilds;
`PDO_REPLOT_ONLY=1` replots without loading the large Seurat input. These flags
cannot both be set. Small sample size, singleton-batch imbalance, response
extreme selection, and potentially correlated clinical features make results
exploratory; no broad clinical or causal conclusion is warranted.

References: the [Bioconductor multi-sample single-cell workflow](https://www.bioconductor.org/books/3.15/OSCA.multisample/multi-sample-comparisons.html),
the [edgeR user guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf),
and the [CAMERA method](https://doi.org/10.1093/nar/gks461).
