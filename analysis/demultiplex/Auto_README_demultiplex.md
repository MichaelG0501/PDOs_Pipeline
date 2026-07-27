# Auto PDO Demultiplex Rerun

This folder contains the organized rerun and audit scripts for the multiplexed PDO pools.

Inputs are read from:

- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/X204SC25083484-Z01-F001/X204SC25083484-Z01-F001/01.RawData`
- `/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/strelka/NT_<donor>/`

Rerun outputs are written to:

- `/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/demultiplex` (for final outputs)
- `/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/demultiplex_intermediate` (for intermediate BAMs and Souporcell output)

Run order:

```bash
qsub -v pool=PDOs_Untreated analysis/demultiplex/Auto_01_cellranger_pdo_pool.sh
qsub -v pool=PDOs_Treated analysis/demultiplex/Auto_01_cellranger_pdo_pool.sh

qsub -v pool=PDOs_Untreated,k=6 analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh
qsub -v pool=PDOs_Treated,k=4 analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh

qsub -v pool=PDOs_Untreated analysis/demultiplex/Auto_03_reference_and_assign.sh
qsub -v pool=PDOs_Treated analysis/demultiplex/Auto_03_reference_and_assign.sh
```

The default donors are `1070,1090,1072,1121,1141,1181` for `PDOs_Untreated` and `1070,1090,1072,1181` for `PDOs_Treated`.
The assignment step uses a parameterized `genotyping_save.R`-style Pearson correlation workflow, not the older `genotype.sh` overlap-count script.

After genotype assignment has completed, export donor-specific count CSVs:

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/demultiplex/Auto_04_write_demultiplexed_counts.R --pool PDOs_Untreated
Rscript analysis/demultiplex/Auto_04_write_demultiplexed_counts.R --pool PDOs_Treated
```

Verify the current pipeline objects and, when the rerun exists, compare old barcode-to-donor assignments against the rerun:

```bash
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/demultiplex/Auto_05_verify_existing_assignments.R
```

The old live scripts referenced `/rds/general/.../PDOs_new`, but that directory is not currently present, so the old `clusters.tsv` and `Genotype_ID_key.txt` for `PDOs_Untreated` cannot be audited directly from their original intermediate files.

To submit the full dependency chain for both pools:

```bash
bash analysis/demultiplex/Auto_00_submit_demultiplex_rerun.sh
```

####################

## Generic single-pool submission and temporary labels

`Auto_00_submit_demultiplex_pool.sh` accepts any folder containing paired
FASTQ files. It retains the historical untreated/treated defaults, but permits
a newly delivered pool to be run without copying or moving FASTQs.

For a pool with WES/VCF references available through the existing donor
configuration, use `reference` mode. For a pool with no WES/VCF yet, use
`temporary` mode: Cell Ranger and Souporcell run normally, only singlet cells
are exported, and each count matrix is named from the observed Souporcell
cluster, e.g. `TEMP_new4samples_SouporcellCluster0_PDO.csv`. Those labels are
not donor identities and must be replaced only after reference genotyping.

The intended new-four-PDO command is:

```bash
bash analysis/demultiplex/Auto_00_submit_demultiplex_pool.sh \
  new4samples \
  /rds/general/project/tumourheterogeneity1/live/ITH_sc/new4samples \
  4 \
  temporary \
  /rds/general/project/tumourheterogeneity1/live/ITH_sc/new4samples/00_counts_matrix_all
```

Heavy Cell Ranger and Souporcell products are kept under
`ephemeral/PDOs_Pipeline/PDOs_outs/demultiplex_intermediate/`. Durable
assignment audits, submission records, and CSV exports are written under live
storage; the explicit final count-matrix directory above is the requested
downstream-compatible export location.
####################
