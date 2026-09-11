# PDO CNV and Numbat methodology

## Scope and interpretation

These are optional CNA diagnostics for malignant PDOs. They do not define the
canonical five PDO states or refine metaprograms. Current CNA scripts may use
the centred state/MP objects only for descriptive overlays; all `legacy_`
state/MP-coupled CNV scripts retain the superseded uncentred route for
provenance and must not feed current analysis.

## Expression-derived InferCNA

`Auto_PDO_infercna.R` estimates expression-derived CNA signal from post-QC PDO
objects using the Carroll 2023 non-malignant epithelial reference and the hg38
Gencode v27 gene order. The main inputs are:

- `PDOs_outs/by_samples/<sample>/<sample>.rds`;
- `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Carroll_2023_reference.rds`;
- `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`.

`SUR843T3_PDO` remains excluded. Heatmap colour clipping uses the 98.5th
percentile of the absolute non-zero CNA values so a few extreme cells do not
flatten the slide-scale display. Reference-derived diagnostic cutoffs for CNA
correlation and signal are the reference mean plus two standard deviations.
Those cutoffs are QC guides, not validated malignant/benign classifiers; all
PDO cells are already malignant by experimental design.

`Auto_PDO_cnv_compare.R` compares available expression-derived and Numbat
profiles. Expression values are rescaled for visual comparability and smoothed
with a rolling window of at most 100 ordered genes. Genome-wide agreement is
also summarized in fixed bins. The default gain/loss classification tolerance
is 0.05 around neutrality. These are visualization/technical concordance
choices and must not be interpreted as clinical copy-number calls.

## Numbat inputs and primary model

`Auto_PDO_numbat_export_inputs.R` writes the persistent live manifest, raw
count matrices, cell maps, and barcode/BAM audit information under
`PDOs_outs/Auto_PDO_numbat/`. The master wrapper and all model wrappers execute
from the live repository. The container image, local R library, source copy,
Singularity cache, and temporary files are reproducible heavy caches under:

`/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_PDO_numbat/`

`Auto_run_pdo_numbat_pileup.sh` uses the official container's
`pileup_and_phase.R` with hg38, the bundled 1000 Genomes panel, automatic UMI
tag detection, and the `CB` cell tag. `Auto_PDO_numbat_run_sample.R` requires at
least 50 count/allele-overlapping cells and runs `run_numbat` with:

- `max_iter = 2`;
- `t = 1e-5`;
- `gamma = 20`;
- `init_k = 3`;
- `min_cells = 50`;
- clonal LOH calling and convergence checks enabled.

These are explicit analysis parameters, not data-adaptive significance
thresholds. Changing them creates a different clone model and requires a
forced rebuild plus a methodology update. Raw Numbat model objects may be
large, but manifests, barcode mappings, clone/joint posterior summaries,
consensus segments, completion records, and any input needed by another
script must remain in live storage.

## Conservative tree re-cut

`Auto_PDO_numbat_conservative_recut.R` is an optional sensitivity layer that
does not overwrite primary Numbat output. By default it considers tree cuts up
to three (`PDO_NUMBAT_CONSERVATIVE_N_CUT=3`) and retains a clone only when it
contains at least the larger of 20 cells or 3% of the sample. Smaller clones
are reassigned to the retained clone with the greatest posterior probability;
if those posterior columns are unavailable, the largest retained clone is the
documented fallback. The script exports the full cut sweep so the selected cut
and every minimum-size decision can be audited.

The 20-cell floor prevents very small posterior groups from being presented as
stable clones; the 3% criterion scales that guard to sample size. These are
conservative visualization/sensitivity choices, not universal biological
truths. Override variables must be recorded in the run summary and figure
caption.

## Outputs, limitations, and reruns

Current outputs use live `intermediate/`, `tables/`, `figures/`, `logs/`, or
workflow-specific `by_samples/` folders. InferCNA is affected by transcription,
cell cycle, gene coverage, and reference choice. Numbat depends on allelic
coverage, phasing, and model convergence. Agreement between them strengthens a
technical observation but does not replace DNA-based validation.

All heavy runs use PBS with `#PBS -koed`. Cached model reuse is allowed only
when the input manifest, parameters, and algorithm version match; use
`PDO_FORCE_REBUILD=1` where supported after changing them.
