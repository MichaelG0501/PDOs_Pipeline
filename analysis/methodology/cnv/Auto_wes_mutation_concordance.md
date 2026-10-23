# Existing-call OAC PDO–tumour mutation plots

Script: `analysis/cnv/Auto_PDO_wes_mutation_concordance.R`.
This terminal workflow adapts panels 2a and 2b of Li et al., *Nature
Communications* 9, 2983 (2018), https://doi.org/10.1038/s41467-018-05190-9.
OAC and EAC describe the same adenocarcinoma histology here.

## Samples and inputs

The source is
`/rds/general/project/spatialtranscriptomics/live/sarek_mutect/`.
Use Mutect2 outputs for PDOs 1070, 1072, 1090, 1121, 1141 and 1181,
and primary tumours 1072, 1090 and 1121. Each was called against its
patient-matched NT sample. PDOs 629 and 727 have no matched normal in this
run and are excluded. The normal is the calling control, not a plotted
tumour sample. Read raw filtered VCFs under `variant_calling/mutect2/`
and corresponding VEP VCFs under `annotation/mutect2/`.

## Selection and annotation

Retain only records with `FILTER=PASS`. Do not add depth, VAF, population
frequency or caller-consensus thresholds. Split multiallelic records into
alternate alleles and retain the corresponding Mutect2 FORMAT AF, AD and DP.
AF is the caller's allele-fraction estimate; it is not the ANNOVAR population
AF column. Match sample columns by their explicit `SUR<patient>_<sample>`
names.

The nonsynonymous coding subset uses VEP protein-coding annotations with
frameshift_variant, stop_gained, stop_lost, start_lost, inframe_insertion,
inframe_deletion, missense_variant or protein_altering_variant. Synonymous,
noncoding and splice-only annotations are excluded from that subset.
If several transcripts annotate one allele in one gene, retain one
consequence in the priority order above. This removes transcript duplication;
it does not choose one variant per gene. Every distinct variant is retained.
VEP's trimmed Allele field links annotations to their alternate allele.
Ambiguous links or missing critical fields stop the script.

Panel a uses the 17 genes displayed in the published figure, in its original
order. Keep genes with no detected mutations as grey rows. Coloured bands show
individual variants, with their AF values to two decimal places. Top bars count
qualifying variant–gene events in those 17 genes. The current-cohort frequency
bar counts affected PDOs out of six; it does not double-count matched tumours
as additional patients. Separate PDO and tumour counts are saved in the
frequency table. A gene's inclusion in this published list does not establish
that each observed variant is a functional driver.

Display the three complete pairs first (1072, 1090, 1121), with tumour then
PDO adjacent within each pair and larger gaps between patients. Display
unpaired PDOs 1070, 1141 and 1181 afterwards. Patient labels span paired
columns rather than repeating the identifier under each sample.

## Published external frequencies

The plotting script hard-codes the fixed gene list,
published MutSigCV asterisks, and approximate Secrier & Li (n=129) and TCGA
OAC (n=89) frequencies. These are descriptive historical references.
They are not newly calculated estimates or formal comparisons against this
small cohort.

The Figure 2 image was downloaded from
https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41467-018-05190-9/MediaObjects/41467_2018_5190_Fig2_HTML.png
on 2026-10-20. Its 2053 × 887 pixel copy is retained at
`PDOs_outs/Auto_wes_mutation_concordance/reference/Auto_Li2018_figure2.png`.
The original 2018 Supplementary Data 2 provides variant allele frequencies,
not the external cohort-frequency table. The 2016 Supplementary Table 4
provides MutSigCV statistics, not all frequencies required here.

With user authorisation, digitise each horizontal black fraction against
its grey background. Secrier bars span x=918–999; TCGA bars span x=1069–1150
(zero-based image pixels). Row centres are y=164+37i for i=0,...,16.
For each pixel calculate `clamp((grey - mean(R,G,B))/grey, 0, 1)`, sum across
the 82-pixel bar, divide by 82, and multiply by 100. Round to the nearest
percentage point. Pixel resolution is approximately 1.22 percentage points;
edge placement and anti-aliasing add uncertainty. All figures and tables
label these estimates approximate. Do not infer exact mutation counts from
these values. The hard-coded reference table records the source and pixel coordinates per gene.

## Concordance

For each of the three matched pairs, compare exact CHROM:POS:REF:ALT keys.
Report PDO-only, shared, and tumour-only counts divided by the union of called
variants in the pair. Produce a total row using original Mutect2 PASS alleles
and a nonsynonymous row using the VEP coding subset across all genes.
Deduplicate alleles across gene annotations before computing overlap.

This is descriptive callset overlap. There is no pooled recall, forced
genotyping, left alignment, coverage rescue or claim of biological absence.
Equivalent indels represented differently can be counted separately.
Low tumour purity and uneven coverage can affect overlap. The VEP run used
`--filter_common`; therefore the nonsynonymous subset reflects annotations
available in those existing files. Total counts use the original VCFs.
The original Sarek metadata has `wes=false`, `WES=true`, and whole-genome
calling intervals. This workflow does not repair or validate that upstream
configuration. No signatures or copy-number analyses are performed.

## Storage, execution and validation

Run in dmtcp on PBS with one CPU, 8 GB RAM and a 30-minute walltime, through
`job_scripts/cnv/Auto_PDO_wes_mutation_concordance.sh`.
The stable live output root is `PDOs_outs/Auto_wes_mutation_concordance/`.
Reruns regenerate the documented current tables, figures, results index and
R session log in that directory. Retained obsolete attempts belong outside
the active output tree in the ignored project-local session archive. Save all PASS alleles, coding annotations, driver
events, overlap membership, aggregate counts and reference frequencies as
CSV source tables. Export one combined final figure in PDF and PNG.
Save `reports/Auto_results.md` and `logs/Auto_sessionInfo.txt`. The script checks allele
mapping, duplicate keys, VAF bounds, nonempty inputs and proportions summing
to one. Verify PBS completion, source-table arithmetic and the rendered
figure before handoff.
