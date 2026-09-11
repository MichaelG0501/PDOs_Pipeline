# New-four-PDO low-pass WGS Strelka methodology

## Status and purpose

This active upstream workflow produces donor-specific reference genotype VCFs
for the `new4samples` Souporcell pool. It follows the established
nf-core/Sarek 3.4.4 germline Strelka route used by the historical PDO WES
workflow, while limiting variant calling to Strelka.

## Inputs

- Sarek samplesheet:
  `analysis/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_samplesheet.csv`
- Two paired FASTQ lanes for each of `SUR1346`, `SUR1363`, `SUR1384`, and
  `SUR1391`, under:
  `/rds/general/project/spatialtranscriptomics/live/WES_PDO/new4_WGS_X204SC26071494-Z01-F001/X204SC26071494-Z01-F001/01.RawData/`

The samplesheet uses one unique patient and sample ID per donor and repeats the
sample ID across its two lanes so Sarek merges the lanes. `status=0` is a
technical instruction to Sarek to run single-sample germline Strelka; it does
not redefine the biological PDO material as normal tissue.

## Method

The PBS wrapper runs nf-core/Sarek `3.4.4` with:

- `-profile singularity`
- `--genome GATK.GRCh38`
- `--trim_fastq`
- `--tools strelka`
- no `--wes` flag, because the delivery is low-pass whole-genome sequencing

Sarek performs its normal mapping and preprocessing route, merges the two lanes
per donor, and generates single-sample Strelka germline output. The downstream
demultiplex workflow requires the variant-only VCF:
`<sample>.strelka.variants.vcf.gz` and its tabix index.

The documented Souporcell/Assign_Indiv_by_Geno Pearson genotype-correlation
method is not changed by this workflow.

## Outputs and storage

Heavy Nextflow work files and the full Sarek result tree remain under:

`/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/demultiplex_intermediate/strelka_sarek/new4samples/`

The required persistent VCFs and indexes are copied without overwriting to:

`/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/strelka/SUR<donor>/SUR<donor>.strelka.variants.vcf.gz`

A lightweight run manifest is written to:

`PDOs_outs/demultiplex/strelka/logs/Auto_new4_lowpass_wgs_sarek_strelka_run_summary.tsv`

These VCFs are downstream inputs, not terminal outputs. Their shared-site
counts and correlation matrix must be reviewed after genotype assignment
because the WGS delivery is approximately 1x coverage.

## Submission

```bash
qsub analysis/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_strelka.sh
```
