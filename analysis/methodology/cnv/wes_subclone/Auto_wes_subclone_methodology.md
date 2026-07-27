# Auto_wes_subclone_methodology.md
####################

# WES FACETS -> PyClone-VI Methodology

This workflow builds bulk tumour-normal subclone inputs for the two PDO samples
with enough PASS Mutect2 SNVs for higher-confidence inference:

- `PDO_1090_vs_NT_1090`
- `PDO_1181_vs_NT_1181`

The default submitter deliberately excludes:

- `PDO_1070_vs_NT_1070`
- `PDO_1072_vs_NT_1072`
- `PDO_1121_vs_NT_1121`
- `PDO_1141_vs_NT_1141`

The excluded samples can be run only by manually setting
`ALLOW_LOW_CONFIDENCE=1`; this is for sensitivity checks, not the default
analysis.

## Software Environment

`Auto_prepare_wes_subclone_env.sh` creates a local mamba environment at
`PDOs_outs/Auto_wes_subclone/conda_env` with `r-facets`, `pyclone-vi`,
`samtools`, and `bcftools`. It also verifies that `snp-pileup` is available,
compiling the source bundled with the FACETS R package against the local htslib
if Bioconda does not put an executable directly on `PATH`.

## Reference Guard

`Auto_validate_sarek_reference.sh` parses the original Sarek launch script:

`/rds/general/project/spatialtranscriptomics/live/WES_PDO/sarek/sarek.sh`

The launch used nf-core/Sarek `3.4.4`, Singularity, `--genome GATK.GRCh38`,
`--WES`, and `--tools mutect2,strelka,cnvkit,vep`. No explicit `--fasta`,
`--igenomes_base`, or `--igenomes_ignore` was passed in the launch command.

Before FACETS starts, the candidate FASTA from `WES_SUBCLONE_REF_FASTA` must:

- exist and have an existing `.fai`;
- match tumour and normal CRAM `@SQ` contig names (`SN`);
- match contig lengths (`LN`);
- match MD5 values (`M5`) wherever both CRAM and FASTA dictionary provide M5;
- decode a small CRAM slice through `samtools view -T`.

If any check fails, the wrapper stops before `snp-pileup`.

## FACETS

FACETS is run from matched tumour-normal CRAMs using `snp-pileup` and a required
common-SNP VCF supplied as `WES_SUBCLONE_FACETS_SNP_VCF`. `Auto_prepare_facets_resources.sh`
downloads the Broad/GATK `Homo_sapiens_assembly38` reference and converts UCSC
`snp151Common` into a biallelic, bgzipped, tabix-indexed SNP VCF under
`PDOs_outs/Auto_wes_subclone/resources/`. It also populates an htslib
`REF_CACHE` from the same FASTA so `snp-pileup` can decode CRAMs whose embedded
Sarek work-directory FASTA path no longer exists. The Sarek Mutect2 VCF is not
used as a substitute for FACETS SNP pileups. CNVkit total copy number is not
used as a substitute for allele-specific copy number.

FACETS outputs purity, ploidy, total copy number, minor copy number, and derived
major copy number. PDO samples are patient-derived organoids and are expected to
be mostly tumour, but tumour content is still taken from FACETS rather than set
manually to 1.

## PyClone-VI

PyClone-VI input is built only from Mutect2 PASS somatic SNVs. Indels are
excluded. `chrX` and `chrY` are excluded by default. For each mutation, the input
uses:

- tumour ref/alt counts from Mutect2 `AD`;
- `major_cn`, `minor_cn`, and `normal_cn` from FACETS segments;
- `tumour_content` from FACETS purity.

The PyClone-VI wrapper writes per-sample inputs, mutation maps, model output,
and cluster results under `PDOs_outs/Auto_wes_subclone/`.

## Visualisation And Reliability Assessment

`legacy_Auto_plot_wes_subclone_results.R` is the terminal QC/interpretation script. It
uses PyClone-VI CCF and assignment-probability outputs together with FACETS
purity, ploidy, and allele-specific segments. It does not call clone number from
raw VAF alone.

The summary PDF contains one page per WES pair:

- PyClone-VI mutation-level CCF by cluster;
- raw tumour VAF versus copy-number-corrected CCF;
- mutation assignment probability by cluster, with 0.8 marked as a high-confidence
  reference threshold;
- mutation support per cluster;
- FACETS autosomal total copy number segments, highlighting minor-copy-zero/LOH
  segments;
- a short reliability text panel.

The reliability table records `n_pyclone_variants`, `n_pyclone_clusters`,
cluster sizes, assignment probabilities, minimum median-CCF separation, and
FACETS purity/ploidy. The current conservative calls are:

- `strong` if variant count, smallest cluster size, assignment probability, and
  median-CCF separation all pass stringent thresholds;
- `moderate` if support is adequate but one or more stringent thresholds are not
  met;
- `weak/provisional` if the exact cluster count should not be treated as
  definitive without sensitivity/re-run support.

The script also writes an optional scRNA Numbat context table when
`PDOs_outs/Auto_PDO_numbat_subclone_mp_conservative/Auto_PDO_numbat_subclone_summary.csv`
is present. This is used as biological context only; WES clone reliability is
called from the WES/PyClone-VI diagnostics.

`legacy_Auto_pyclone_sensitivity.R` reruns PyClone-VI on the same input tables with
maximum-cluster caps of 5, 10, 20, and 40 using three restarts per cap. This is
a stability screen rather than a replacement for the main 40-cluster/10-restart
run. In the current run:

- `PDO_1090_vs_NT_1090` remains at 2 occupied clusters across all tested caps;
- `PDO_1181_vs_NT_1181` gives 2 occupied clusters at cap 5 but 3 occupied
  clusters at caps 10, 20, and 40, with low assignment probabilities in the
  split higher-CCF clusters.

Therefore the current WES evidence supports a more stable two-cluster solution
for `PDO_1090_vs_NT_1090`, while the exact three-cluster solution for
`PDO_1181_vs_NT_1181` should remain provisional.

## High-Resolution WES/scRNA CNA Audit

`legacy_Auto_wes_scrna_subclone_highres_audit.R` was added after reviewing the first
WES/scRNA matching figures. It diagnoses two technical issues:

- the original generated WES bulk/subclone `.cns` files were FACETS-resolution
  (36 segments for `PDO_1090_vs_NT_1090`, 32 for `PDO_1181_vs_NT_1181`) even
  though Sarek CNVkit contains higher-resolution `.cns`/`.somatic.call.cns`
  segment grids;
- concordance calculations should bin segmented CNA profiles by segment overlap
  at the bin midpoint, not by testing whether segment boundary points fall
  inside each bin.

The high-resolution audit writes new outputs under live
`PDOs_outs/Auto_wes_subclone/tables/cns_highres/`,
`figures_highres/`, `tables/visualisation_highres/`, and
`PDOs_outs/cnv/cnv_compare_highres/`. It reads large Numbat/inferCNA inputs
from the ephemeral project copy but does not copy those intermediate files into
live.

The WES per-subclone CNA profiles remain projections from bulk CNV data:
PyClone-VI clusters SNVs, while FACETS/CNVkit segment CNAs in bulk. Without a
joint phylogenetic model of SNVs and CNAs, WES alone cannot assign every CNA
breakpoint uniquely to a PyClone SNV cluster. High-resolution projected profiles
are therefore useful for checking whether a WES SNV cluster is compatible with
scRNA CNA clones, but they should not be interpreted as independently inferred
single-subclone CNA genomes.

Numbat clone profiles must be plotted on their native `log2(phi_mle_roll)`
scale when interpreting copy-number state, because `phi_mle_roll` is the total
copy-number ratio relative to diploid. Median-centering a clone can make WES and
Numbat shapes look more similar, especially for polyploid samples such as
`SUR1181_Treated_PDO`, but that centered view is only a shape/breakpoint
comparison and erases global aneuploidy. Label any centered output explicitly as
baseline-centered/shape-only; do not call it the native per-clone CNA profile.

## Absolute Bulk WES CNA Scaling

`legacy_Auto_wes_absolute_cna_compare.R` writes a scale-corrected WES bulk CNA track
from FACETS allele-specific segments. The exported value is
`log2(total_cn / 2)`, i.e. absolute tumour total copy number scaled to a diploid
baseline. This is intentionally different from the original Sarek CNVkit `.cns`
log2 ratio, which is centered for relative depth visualization and can make a
polyploid sample appear visually close to zero.

For the current run, `PDO_1181_vs_NT_1181` has FACETS purity 0.663, ploidy
3.206, median total copy number 3, and median diploid-scaled log2 0.585. The
absolute WES track is therefore broadly amplified, matching the native Numbat
interpretation for `SUR1181` better than the centered CNVkit track. The script
also projects the FACETS absolute values onto the CNVkit segment grid for a
high-resolution `.cns` visualization file, but the absolute copy states still
come from FACETS total copy number.

Native Numbat comparisons in this script use `log2(phi_mle)` or
`log2(phi_mle_roll)` directly. They are not median-centered. Correlations are
reported at 1 Mb and 5 Mb bins as scale checks, not as proof that WES PyClone
clusters and scRNA Numbat clones represent the same biological clone.

## THetA2 Clone-Specific WES CNA Attempt

True clone-specific WES CNA requires a model that deconvolves copy-number
segments into tumour subpopulations. PyClone-VI cannot do this from its output:
it consumes major/minor copy number and tumour content as inputs for SNV
clustering, then outputs SNV cluster cellular prevalence and assignment
probabilities.

The implemented model-based route is:

1. `Auto_prepare_theta2_clone_cna_env.sh` builds separate ephemeral conda
   environments under `ephemeral/PDOs_Pipeline/PDOs_outs/Auto_wes_clone_cna/`:
   `theta2_env` for the Python 2.7-era THetA2 package and `cnvkit_py310_env`
   for a Python 3 CNVkit install pinned to Python 3.10/pandas <2. A combined
   env is not used because `theta2` forces Python 2.7 while current CNVkit
   imports Python 3-only symbols.
2. `Auto_run_theta2_clone_cna_sample.sh` uses the conditional-shift WES bulk
   `.cns` when present, falling back to the original Sarek CNVkit tumour `.cns`,
   plus `reference.cnn` to run `cnvkit.py export theta`.
3. `Auto_make_theta2_snp_counts.R` converts the existing FACETS `snp-pileup`
   into THetA2 tumour/normal SNP count files by retaining autosomal common SNPs
   with normal depth >= 20, tumour depth >= 20, and matched-normal allele
   fraction 0.25-0.75.
4. `RunTHetA` is run with `--BAF`, and `cnvkit.py import-theta` imports the
   THetA2 `BEST.results` file back into CNVkit `.cns` format. Imports are staged
   in a per-job ephemeral directory and only those newly imported files are
   copied into the live output directory. This avoids cross-sample `.cns` files
   being re-prefixed and mixed into later imports.

This workflow is the first proper WES clone-CNA attempt in this repository. It
may return fewer segments than the input CNVkit `.cns` because THetA2 performs
model selection/significance filtering on segments; that reduction is expected
for a true clone-CNA model and should not be confused with the earlier
resolution loss from PyClone projection. The current high-confidence runs
completed for `PDO_1090_vs_NT_1090` and `PDO_1181_vs_NT_1181`, selected `n=2`
BEST solutions, and imported one tumor `.cns` per sample (159 and 104 regions,
respectively). Forced `n=3` imports from the corrected conditional-shift input
produced two clone `.cns` files per sample. SUR1090 clone 2 was moderately
concordant with Numbat (5 Mb correlations up to about 0.61), but neither
THetA2 clone outperformed the conditional-shift CNVkit bulk track. SUR1181
THetA2 clone 2 preserved more amplification than HATCHet but remained coarse
and only partly concordant with Numbat. These outputs are therefore cautious
model-based sensitivity checks rather than definitive high-resolution clone
genomes.

####################
## 2026-07-08 High-Resolution Absolute CNVkit Scaling

`Auto_wes_absolute_highres_subclone_compare.R` supersedes the first FACETS-only
absolute comparison for presentation figures. The first absolute comparison was
biologically scaled but visually too coarse because it plotted FACETS segment
breakpoints directly. The corrected script reads the CNVkit-resolution WES
profiles from `PDOs_outs/Auto_wes_subclone/tables/cns_highres/`, preserves the
original CNVkit segment grid, and adds the FACETS ploidy offset:

`absolute_log2 = centered_cnvkit_log2 + log2(FACETS_ploidy / 2)`.

This gives the intended absolute WES visualization for polyploid samples while
keeping the high-resolution shape track needed for WES/scRNA matching. For the
current outputs, `PDO_1181_vs_NT_1181` has FACETS ploidy 3.206 and therefore a
log2 offset of 0.681; the corrected bulk CNVkit `.cns` has median absolute log2
0.664 and 100/106 positive segments. The corrected figures include:

- centered CNVkit WES shape track;
- conditional-shift absolute WES track;
- conditional-shift projected WES SNV-cluster/subclone tracks;
- native Numbat pseudo-bulk and clone CNA tracks;
- inferCNA mean track when available;
- 5 Mb WES/scRNA correlation heatmap.

The current adjusted `.cns` files are written under
`PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_<sample>_*_conditional_shift_absolute.cns`.
The figure filenames intentionally remain
`PDOs_outs/Auto_wes_absolute_cna/figures/Auto_wes_absolute_cna_compare_<sample>.pdf/.png`
so downstream slide links use the corrected high-resolution views.

Important interpretation boundary: these WES subclone rows are still projected
bulk CNVkit profiles stratified by PyClone-VI SNV clusters. They are useful
WES/scRNA compatibility and visualization tracks, but they are not independent
true clone-specific CNA genomes. The THetA2 route remains the model-based WES
clone-CNA attempt; current THetA2 solutions are too low-resolution for the
presentation matching figure.
####################

####################
## 2026-07-09 Integrated PhyloWGS SNV+CNA Route

HATCHet and THetA2 did not produce convincing clone-specific CNA tracks for the
current single-WES samples. The preferred integrated follow-up is therefore:

1. FACETS defines tumour purity, ploidy, and allele-specific CNA segments.
2. Mutect2 PASS SNVs are mapped to those FACETS segments and used as
   copy-number-aware PyClone-VI input.
3. `Auto_make_phylowgs_inputs.R` converts the FACETS/PyClone handoff into
   PhyloWGS `ssm_data.txt` and `cnv_data.txt`, with audit tables preserving
   the PyClone cluster labels and FACETS segment coordinates.
4. `Auto_run_phylowgs_sample.sh` runs PhyloWGS with 4 chains, 1000 burn-in
   samples, and 2500 MCMC samples. The old Python 2 PhyloWGS environment and
   bulky MCMC intermediates are kept under the ephemeral WES subclone tree.
5. `Auto_summarise_phylowgs_results.R` extracts the highest-density posterior
   tree, joins PhyloWGS population assignments back to FACETS CNV events and
   PyClone-mapped SSMs, and writes inherited clone-CNA event tables.

Completed high-confidence runs:

- `PDO_1090_vs_NT_1090`: 154 SSMs, 35 FACETS CNA events, 2500 posterior trees.
  Top tree `712` has 6 populations including root, 5 event populations, 35
  assigned CNV events, and 154 assigned SSMs. Population 1 is the trunk-like
  event population (cellular prevalence 0.782, 112 SSMs, 26 CNVs), with
  descendant populations carrying additional CNA/SNV events.
- `PDO_1181_vs_NT_1181`: 99 SSMs, 31 FACETS CNA events, 2500 posterior trees.
  Top tree `944` has 6 populations including root, 5 event populations, 31
  assigned CNV events, and 99 assigned SSMs. Population 1 is trunk-like
  (cellular prevalence 0.659, 25 SSMs, 12 CNVs), followed by CNA-bearing
  descendant populations.

Final live outputs:

- `PDOs_outs/Auto_wes_subclone/reports/phylowgs/<sample>/<sample>.trees.zip`
- `PDOs_outs/Auto_wes_subclone/reports/phylowgs/<sample>/<sample>.summ.json.gz`
- `PDOs_outs/Auto_wes_subclone/reports/phylowgs/<sample>/<sample>.muts.json.gz`
- `PDOs_outs/Auto_wes_subclone/reports/phylowgs/<sample>/<sample>.mutass.zip`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs/Auto_phylowgs_top_tree_summary.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_phylowgs_population_top_tree.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_phylowgs_cnv_assignment_top_tree.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_phylowgs_ssm_assignment_top_tree.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/Auto_<sample>_phylowgs_clone_cna_inherited_top_tree.csv`

Interpretation boundary: PhyloWGS is the correct model class for integrating
SNVs and CNAs into one phylogeny, and it is a better match to the stated goal
than PyClone-only projections, HATCHet, or THetA2 for these samples. However,
single-WES data still cannot guarantee a unique full CNA genome for every
SNV-defined clone. The inherited clone-CNA tables should be interpreted as the
highest-density PhyloWGS assignment of FACETS CNA events to tree populations,
not as independently validated complete clone genomes.
####################

####################
## 2026-07-08 Conditional Baseline Shift Update

The high-resolution absolute comparison was updated after review of the
SUR1090 figures. The FACETS ploidy offset is no longer applied blindly to every
WES pair. The script now calculates the native Numbat median log2 baseline from
gene-level `bulk_clones` `phi_mle_roll` values across the matched scRNA samples:

- if the native Numbat median is below 0.25 log2, the applied WES shift is 0 and
  the CNVkit-centered shape is retained for the absolute/display row;
- if the native Numbat median is at least 0.25 log2, the applied WES shift is
  `log2(FACETS ploidy / 2)`.

Current result: `PDO_1090_vs_NT_1090` has native Numbat median -0.036, so the
applied shift is 0 despite a FACETS ploidy offset of 0.919. `PDO_1181_vs_NT_1181`
has native Numbat median 1.057, so the applied shift is 0.681. This keeps SUR1090
from being artificially shown as globally amplified while preserving the SUR1181
polyploid amplification display.

Numbat pseudo-bulk is also plotted from a gene-level weighted `bulk_clones`
profile when available, not from the lower-resolution consensus-segment table.
The current adjusted WES `.cns` files use the suffix
`*_conditional_shift_absolute.cns`.
####################

####################
## 2026-07-08 HATCHet Clone-CNA Audit

HATCHet was added as the preferred model-based WES clone-CNA audit after the
PyClone-VI projection and THetA2 outputs failed to provide convincing
clone-specific CNA profiles. The implemented route is:

1. `Auto_prepare_hatchet_clone_cna_env.sh` installs Bioconda HATCHet under the
   ephemeral WES clone-CNA output tree.
2. `Auto_prepare_hatchet_cbc_solver.sh` installs CBC into the same environment.
   The packaged C++ solver attempts to link a missing Gurobi runtime on CX3, so
   `Auto_run_hatchet_clone_cna_sample.sh` writes a local `hatchet.ini` forcing
   Pyomo/CBC.
3. `Auto_make_hatchet_bb.R` builds 250 kb HATCHet `.bb` inputs from Sarek
   CNVkit `.cnr` read-depth ratios and FACETS germline-heterozygous SNP pileups.
   Read depth is scaled by the same conditional shift used for WES visualization
   (`2 ^ (cnr_log2 + applied_shift_log2)`), so SUR1181 no longer enters HATCHet
   as a centered near-diploid profile.
4. `Auto_run_hatchet_clone_cna_sample.sh` runs `cluster-bins`, then
   `compute-cn` with a bounded clone range. The current default is `2,4`, and
   FACETS purity is passed through HATCHet's purity option (`-P`) when available.
5. `Auto_hatchet_clone_cna_compare.R` parses HATCHet `best.bbc.ucn`, converts
   allele-specific `A|B` calls into `log2((A+B)/2)`, and compares both the pure
   clone genotype and the mixture-expected bulk profile to conditional-shift WES
   bulk and native high-resolution Numbat `bulk_clones`.

Current HATCHet selected solutions:

- `PDO_1090_vs_NT_1090`: diploid `n=2`, selected `best.bbc.ucn` contains one
  tumour clone column (`u_clone1=0.801`, `u_normal=0.199`), median clone log2
  0.000, mean clone log2 -0.161. The corrected run used FACETS purity 0.801 and
  no conditional read-depth shift.
- `PDO_1181_vs_NT_1181`: diploid `n=2`, selected `best.bbc.ucn` contains one
  tumour clone column (`u_clone1=0.663`, `u_normal=0.337`), median clone log2
  -1.000, mean clone log2 -0.693. The corrected run used FACETS purity 0.663
  and conditional read-depth shift 0.68085. HATCHet did explore tetraploid
  candidates after the shift but still chose the diploid model by its internal
  model-selection criterion.

Interpretation: HATCHet did infer allele-specific absolute CNA and a tumour
clone proportion, but model selection did not support multiple tumour CNA
clones in either single-WES sample. The current single HATCHet clone is not a
plotting omission: the selected `best.bbc.ucn` files contain only
`cn_clone1/u_clone1`. For SUR1090 the conditional-shift CNVkit bulk remains more
concordant with Numbat than the HATCHet clone profile. For SUR1181 HATCHet's
selected clone is biologically discordant with the globally amplified Numbat and
conditional-shift CNVkit profiles; it should be treated as a failed/low-support
clone-CNA deconvolution for this sample. The conditional-shift CNVkit bulk
remains the clearest whole-WES visualization.

Cleanup: obsolete PyClone projection figures/tables, low-resolution `.cns`,
old globally ploidy-adjusted `.cns`, copied HATCHet intermediate files, and
root PBS debug logs were removed. The exact removal manifest is
`PDOs_outs/Auto_wes_clone_cna/tables/Auto_wes_cleanup_manifest.tsv`.
####################

####################
## 2026-07-09 PhyloWGS/Numbat Visualization

`Auto_plot_phylowgs_numbat_compare.R` was added to plot the integrated
PhyloWGS result directly against native high-resolution Numbat clone profiles.
It reads the top-tree inherited clone-CNA event tables from
`PDOs_outs/Auto_wes_subclone/tables/phylowgs/<sample>/`, converts FACETS total
copy number to `log2(total_cn / 2)`, and plots each PhyloWGS event population
against Numbat `bulk_clones_final.tsv.gz` gene-level `log2(phi_mle_roll)`.

Final visualization outputs are:

- `PDOs_outs/Auto_wes_subclone/figures/phylowgs/Auto_phylowgs_clone_cna_compare_<sample>.pdf/.png`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_visualisation_summary.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_correlations.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_segments_for_plot.csv`

The plots are intentionally diagnostic. PhyloWGS assigns FACETS CNA events to
phylogenetic populations; it does not create a dense genome-wide CNA track like
Numbat. Missing WES event bins are therefore treated as neutral only for the
correlation calculation. In the completed SUR1090/SUR1181 figures, the
PhyloWGS event tracks are sparse and poorly correlated with Numbat, while the
conditional-shift CNVkit bulk row remains the more faithful whole-genome WES
CNA visualization.
####################

####################
## 2026-07-10 PhyloWGS Genome-Order Correction And Final Clone Plotting

The 2026-07-09 sparse PhyloWGS/Numbat visualization was superseded after
auditing the FACETS pileups and plot semantics. The root input problem was
that the common SNP VCF used by `snp-pileup` was lexicographically sorted
(`chr1`, `chr10`, ...) rather than in reference contig order. That truncated
the effective FACETS/PyClone evidence to the first ordered chromosomes in the
CRAM traversal and led to chromosome-sparse PhyloWGS CNA event plots.

`Auto_prepare_facets_snp_genome_order.sh` now creates the corrected resource:

- `PDOs_outs/Auto_wes_subclone/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.genome_order.vcf.gz`

FACETS and PyClone-VI were force-rerun for both high-confidence WES pairs using
this genome-order VCF and the ephemeral reference/env paths. The corrected
FACETS pileups are genome-wide:

- `PDO_1090_vs_NT_1090`: 1,030,713 SNP rows, FACETS ploidy 3.409,
  PyClone-VI 1,096 variants and 3 clusters.
- `PDO_1181_vs_NT_1181`: 1,027,554 SNP rows, FACETS ploidy 3.220,
  PyClone-VI 689 variants and 2 clusters.

PhyloWGS was rerun from fresh MCMC directories using
`PHYLOWGS_RUN_SUFFIX=_genome_order_20260709`, then summarized again. The new
top trees assign CNA events across all autosomes:

- `PDO_1090_vs_NT_1090`: 2,500 posterior trees, top tree 137, 6 event
  populations, 284 assigned CNVs, 1,096 assigned SSMs.
- `PDO_1181_vs_NT_1181`: 2,500 posterior trees, top tree 2353, 5 event
  populations, 223 assigned CNVs, 689 assigned SSMs.

The terminal visualization script was also corrected. The PhyloWGS population
rows no longer copy the WES CNVkit bulk backbone. They now show the final
PhyloWGS clone CNA profile: inherited FACETS/PhyloWGS CNA events are plotted on
the conditional WES baseline, and non-event intervals are filled as neutral.
The top WES CNVkit row remains only as a reference. Baseline alignment uses
`log2(total_cn / 2) - facets_ploidy_shift_log2 + applied_shift_log2`, so
SUR1090 is not artificially pushed into a globally amplified state, while
SUR1181 keeps the amplified absolute scale.

The Numbat loader was also reverted to prefer native high-resolution
`bulk_clones_final.tsv.gz` from the ephemeral sample output. Conservative live
Numbat files are only a fallback and are treated as already log2-scaled when
used, avoiding double-log transformation artifacts.

The refreshed terminal outputs are:

- `PDOs_outs/Auto_wes_subclone/figures/phylowgs/Auto_phylowgs_clone_cna_compare_<sample>.pdf/.png`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_visualisation_summary.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_correlations.csv`
- `PDOs_outs/Auto_wes_subclone/tables/phylowgs_visualisation/Auto_phylowgs_clone_cna_segments_for_plot.csv`
####################
