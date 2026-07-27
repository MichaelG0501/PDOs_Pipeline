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
