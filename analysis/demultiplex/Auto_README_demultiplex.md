# Auto PDO Demultiplex Rerun

This folder contains the organized rerun and audit scripts for the multiplexed PDO pools.

Inputs are read from:

- `/rds/general/project/tumourheterogeneity1/live/ITH_sc/X204SC25083484-Z01-F001/X204SC25083484-Z01-F001/01.RawData`
- `/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/strelka/NT_<donor>/`

Rerun outputs are written to:

- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex`

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
