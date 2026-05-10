#!/bin/bash
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -l walltime=2:00:00
#PBS -N Auto_PDO_CurrentBarcodeGeno
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
cd "$WD"
Rscript analysis/demultiplex/Auto_06_current_barcode_genotype_check.R

echo $(date +%T)
