#!/bin/bash
#PBS -l select=1:ncpus=8:mem=160gb
#PBS -l walltime=24:00:00
#PBS -N Auto_PDO_scVelo
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate velocity

export PYTHONPYCACHEPREFIX="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_velocity_PDO/tmp_pycache"

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
cd "$WD"

python analysis/trajectory/Auto_scvelo_pdo_visualise.py

echo $(date +%T)
