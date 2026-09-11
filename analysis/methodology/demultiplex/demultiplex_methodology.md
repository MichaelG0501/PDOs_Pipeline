# PDO Demultiplex Methodology

This document describes the organized demultiplex rerun workflow under
`analysis/demultiplex/`.

## 1. Aim

The demultiplex workflow reruns CellRanger and Souporcell for multiplexed PDO
pools, builds normal-WES donor reference genotypes, assigns Souporcell clusters
to donors, exports donor-specific count matrices, and verifies current object
assignments.

## 2. Core Inputs

- raw FASTQs under
  `/rds/general/project/tumourheterogeneity1/live/ITH_sc/X204SC25083484-Z01-F001/.../01.RawData/<pool>/`
- CellRanger binary:
  `/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-9.0.1/bin/cellranger`
- transcriptome:
  `/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A`
- Demuxafy/Souporcell container from the live demultiplex/multiplexed folders
- normal Strelka VCFs for donor genotype references
- current PDO RDS objects for verification

## 3. Run Order

1. `Auto_01_cellranger_pdo_pool.sh`
2. `Auto_02_souporcell_pdo_pool.sh`
3. `Auto_03_reference_and_assign.sh`
4. `Auto_04_write_demultiplexed_counts.R`
5. `Auto_05_verify_existing_assignments.R`

`Auto_00_submit_demultiplex_rerun.sh` orchestrates this order with PBS
dependencies.

## 4. Output Standards

Demultiplex outputs intentionally include external staging folders under
`/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/demultiplex_intermediate/` for heavy intermediates and `/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/demultiplex/` for final outputs.
Any new scripts must document external writes clearly in their headers and in
`analysis/ANALYSIS_MAP.md`.

####################

## 5. Optional reference genotyping and temporary cluster exports

`Auto_00_submit_demultiplex_pool.sh` runs one multiplexed FASTQ directory. In
`reference` mode it follows the existing Souporcell-cluster to donor-reference
genotype correlation workflow. In `temporary` mode, intended only where WES or
VCF donor references are not yet available, it intentionally skips genotype
assignment and exports only Souporcell singlets.

Temporary matrices are named
`TEMP_<pool>_SouporcellCluster<id>_PDO.csv`. The cluster ID is not a donor ID,
and the temporary label must be replaced after a WES/VCF-backed rerun. For the
new four-sample delivery, Cell Ranger and Souporcell intermediates remain in
the ephemeral demultiplex tree, while the requested downstream CSV export is
`/rds/general/project/tumourheterogeneity1/live/ITH_sc/new4samples/00_counts_matrix_all/`.
####################

PBS wrappers must include `#PBS -koed` for live logs.

####################

## 6. New-four-PDO Strelka donor references

Before assigning the four `new4samples` Souporcell clusters to `SUR1346`,
`SUR1363`, `SUR1384`, and `SUR1391`, run:

```bash
qsub analysis/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_strelka.sh
```

This Sarek 3.4.4 PBS workflow reads the two low-pass WGS lanes per donor from
`Auto_new4_lowpass_wgs_sarek_samplesheet.csv`, runs only Strelka for variant
calling, keeps heavy Nextflow/Sarek outputs under the PDO ephemeral tree, and
copies each required `*.strelka.variants.vcf.gz` plus `.tbi` to the established
live `spatialtranscriptomics/live/sarek_mutect/variant_calling/strelka/` tree.

Detailed methodology:
`analysis/methodology/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_strelka_methodology.md`.

The downstream Souporcell Pearson genotype-correlation method is unchanged.
Once the four VCFs exist, submit the assignment with
`qsub -v pool=new4samples analysis/demultiplex/Auto_03_reference_and_assign.sh`.
The wrapper maps the `new4samples` pool to donors `SUR1346`, `SUR1363`,
`SUR1384`, and `SUR1391`, retains and merges their biallelic PASS heterozygous
Strelka SNP genotypes, and compares GT dosages at
matching CHROM:POS:REF:ALT sites with the four Souporcell cluster genotypes.
The older normal-WES pools retain the historical `DP>=4 && GQ>=20` filter.
That additional threshold is not applied to the approximately 1x new4samples
references because it leaves only 137 shared sites; Strelka PASS remains the
low-pass quality filter, while the downstream Pearson assignment method and
1,000-shared-site guard remain unchanged.

The live assignment directory is
`PDOs_outs/demultiplex/genotype_assignment/new4samples/`. Its documented
visualization is
`Auto_new4samples_ref_clust_pearson_correlation.png`; the same directory holds
the genotype assignment keys, Pearson correlation matrix, overlap-site matrix,
and assignment summary. The visualization is terminal, while the assignment
key is a downstream input to donor count export.
####################

####################

## 7. Historical NT-reference heatmap regeneration

The original `PDOs_Untreated` and `PDOs_Treated` Cell Ranger BAMs, indexes,
and filtered barcode files remain under the live ITH_sc PDO Cellranger tree,
although their Souporcell and genotype-assignment outputs are absent.
`Auto_02_souporcell_pdo_pool.sh` accepts explicit `input_bam` and
`input_barcodes_gz` PBS variables so these retained inputs can be reused in
place. Souporcell outputs remain in the organized ephemeral demultiplex tree.
The unchanged NT-reference assignment then creates persistent correlation
tables, reciprocal keys, summaries, and comparison heatmaps under the live
PDO demultiplex output tree.

The small downstream-critical `clusters.tsv` is additionally copied to
`PDOs_outs/demultiplex/souporcell_assignments/<pool>/Auto_<pool>_clusters.tsv`
after every successful Souporcell run. The count-export script reads this live
copy. Large BAMs, allele matrices, merged VCFs, and `cluster_genotypes.vcf`
remain ephemeral.
####################
