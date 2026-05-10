#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=8:00:00
#PBS -N Auto_PDO_VerifyDemux
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
cd "$WD"
Rscript analysis/demultiplex/Auto_05_verify_existing_assignments.R

echo $(date +%T)
