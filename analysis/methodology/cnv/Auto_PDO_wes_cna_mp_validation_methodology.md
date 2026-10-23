# Six-PDO bulk WES CNA and scATLAS MP comparison

## Scope and sample identity

This exploratory validation uses the existing CNVkit profiles for SUR1070,
SUR1072, SUR1090, SUR1121, SUR1141 and SUR1181. The user confirmed that WES and
RNA sample names identify the same samples and instructed use of existing CNA
outputs despite the original Sarek WES flag issue. Each donor is paired once
with its `SUR<id>_Untreated_PDO` RNA sample. Treated RNA cultures are excluded;
the six DNA measurements are not duplicated across cells or treatment pairs.

The independent unit is a PDO donor/culture. Six donors do not support a claim
of definitive broad validation. The original Sarek parameters contain
`wes=false` and `WES=true`; Sarek rejected uppercase `--WES`, and its CNVkit
configuration selects WGS mode when `wes=false`. These limitations remain
explicit in figures, manifests and reports. This workflow does not rerun or
correct DNA processing and does not label the six profiles quality-validated.
SUR629 and SUR727 are excluded because their existing profiles are far coarser.

## Inputs and RNA scoring

The executable registry lists all canonical inputs. RNA comes exclusively from
the post-QC `PDOs_outs/PDOs_list_PDOs.rds`. Exact `orig.ident` is checked against
each list element; identity is never inferred from cell barcodes. The exact 17
current centred refined scATLAS gene sets must equal the pinned reference
gene-event results. No PDO-specific MP ID is substituted. Gene sets MP2x,
MP11c and MP18a must remain absent.

UCell uses RNA counts, `maxRank=1500`, chunks of 1000 cells and two workers.
Ranking uses all measured RNA genes in each sample. Signature genes absent
from the measured feature universe are omitted and reported; fewer than five
measured genes in any signature stops execution. The arithmetic mean of cell
scores within each sample defines the one MP observation per donor. scATLAS
retains its pinned subclone activity on the reference scale, which can include
negative values; it is not relabeled as raw UCell. TCGA retains its existing bulk
GSVA scores; score units are not pooled across platforms. Cell scores and
sample summaries remain in live storage. A fingerprint of input size and
mtime gates reuse; `PDO_FORCE_REBUILD=1` explicitly rebuilds the score cache.

## DNA measurements and categories

Read only original `PDO_<id>.somatic.call.cns` files. No RNA-guided shift,
Numbat-derived baseline or projected clone CNA is applied. The primary
predictor is the original continuous segment log2 ratio, not absolute CN or
TCGA GISTIC. Twenty curated genes (ten gain and ten loss) are the exact pinned
reference universe. GRCh38 gene coordinates come from `hg38_gencode_v27.txt`,
retaining the first entry per gene as in the reference. Gene coordinates are
treated as 1-based inclusive and converted to overlap with CNVkit's 0-based,
half-open segments. Overlapping segment means are weighted by overlapping
gene base pairs. Less than 90% segment coverage stops execution. This measures
the containing CNA segment; it does not establish a gene-specific focal event.

The five gene boxplot categories use the original CNVkit estimated integer CN:
0, 1, 2, 3 and at least 4. Labels match deep loss, loss, neutral, gain and
amplified for visual comparison, but are explicitly identified as CNVkit CN
estimates, not GISTIC. A plotted gene must have one containing segment and
essentially complete coverage to assign a category. Ambiguous categories stop
the run. Absolute interpretation remains limited by purity, ploidy, centering
and the original processing mode.

Forty-six arms (chr1–chr22 and chrX, p and q) use the same GRCh38 chromosome
sizes and centromere boundaries as the reference workflow. Segment log2 ratios
are weighted by arm-overlap base pairs. Missing arms stop; missing data are
never filled with neutral values. Arm boxplot gain/loss categories use the
reference TCGA threshold retained in its pinned source cache, currently 0.12.
No threshold is fitted to six PDOs. Arm coverage is exported, including sparse
acrocentric arms. Sex-chromosome comparisons are retained to match reference
figures, without sex, purity or ploidy adjustment; they need particular caution.

## Correlation tests and selection

Compute Spearman rho across exactly six samples. Enumerate all 720 permutations
of MP ranks, including the observed arrangement. The two-sided P is the
proportion with absolute rho at least the observed value (tolerance 1e-12).
Tied ranks are retained, and labeled sample permutations count with their
proper multiplicity. Constant predictors or responses have undefined tests.
BH correction uses the complete family of 340 gene/MP or 782 arm/MP tests,
including untestable hypotheses in the family size. Reference rho, P and FDR
remain pinned, using their original tests; no cross-cohort pooled P is fitted.

The ten gene/MP and ten arm/MP boxplot pairs are fixed by the existing reference
top-ten lists. PDO effects never select or reorder pairs. This prevents choosing
only apparently replicated PDO effects. Full families, selected-pair summaries,
direction agreement and plot sources are exported. Correlation of gene CNA
with MP activity does not prove causality. Gene membership in an MP is flagged.

## Figures and interpretation

Gene rectangles have four pages: gain all trends, gain TCGA-significant pairs,
loss all trends and loss TCGA-significant pairs. Arm rectangles have two pages:
all trends and TCGA-significant pairs. Each has three panels: scATLAS, TCGA and
PDO. The gene predictor in scATLAS is gene RNA expression, whereas PDO and
TCGA use DNA CNA. Colour saturates at rho +/-0.25 to match the reference.
An asterisk means TCGA BH FDR below 0.05 with matching nonzero direction in
the displayed cohort. It does not signify PDO significance. Masks use TCGA
FDR alone and retain opposite effects without an asterisk.

Each boxplot PDF contains ten pages, one fixed pair per page and three cohort
panels. scATLAS gene RNA uses type-7 tertiles with ties kept together; TCGA uses
observed GISTIC categories; PDO uses original estimated CN categories. Arms
use cohort-specific reference thresholds. Group counts, all six PDO points and
donor labels are displayed, including empty categories. The subtitle reports
the continuous correlation and full-family FDR. Grouped Wilcoxon tests are
saved separately when both groups contain at least two observations; they are
not substituted for continuous tests or drawn as unsupported brackets. They
are nominal descriptive tests and must not be treated as corrected evidence.
An additional 20-page PDO scatterplot PDF exposes the continuous observations
and leverage of individual donors without fitting a regression line.

## Storage, execution and verification

All outputs are live under `PDOs_outs/Auto_PDO_wes_cna_mp_validation/`, divided
into intermediate, tables, figures, logs and reports. Plot objects and source
data permit `PDO_REPLOT_ONLY=1` without loading RNA. Replot and force-rebuild
flags are mutually exclusive. Original files remain unchanged. No deletion,
source-data relocation or external upload is performed.

Submit `qsub job_scripts/cnv/Auto_PDO_wes_cna_mp_validation.sh` from the project
root. Use gnmf, two CPUs, 64 GB and two hours. The wrapper uses `#PBS -koed`
and absolute project-root stdout/stderr destinations. The main R script runs
directly on a PBS compute node and contains all analytical logic.

Verification checks exact identities, signature universes, sample means,
source-segment MYC values, every Spearman coefficient via a second formulation,
permutation counts and P-value bounds, complete BH families, unique test rows
and fixed reference-selected pairs. The tracked Python helper uses the pinned
PyMuPDF 1.26.7 installation from the reference workflow. It requires five PDFs
with 46 total pages, nonempty text and no text outside page bounds, and writes
five first-page previews. Numerical and PDF reports and job completion must
pass before adopting outputs; representative previews are visually inspected.

## Verified production result, 19 October 2026

PBS job `4412928.pbs-7` completed with exit status zero. The final five PDFs
contain 46 pages; all page-count and text-bound checks passed. First-page
previews were visually reviewed, including the corrected cohort-specific
activity labels. The analysis contains 24,534 RNA cells from six independent
untreated PDO donors. No PDO association passes BH FDR below 0.05 across the
340 gene/MP or 782 arm/MP hypotheses. Five of ten fixed gene pairs and five of
ten fixed arm pairs match TCGA's nonzero effect direction. Direction matching
alone is not statistical replication. For MYC–MP2+, PDO rho is 0.0857 and the
exact two-sided permutation P is 0.9194. These results retain the documented
existing-CNA processing limitations and should be interpreted as exploratory.
