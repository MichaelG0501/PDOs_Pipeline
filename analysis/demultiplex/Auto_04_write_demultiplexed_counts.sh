#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=8:00:00
#PBS -N Auto_PDO_WriteDemuxCounts
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

pool="${pool:-}"
if [[ -z "$pool" ]]; then
  echo "ERROR: submit with -v pool=PDOs_Untreated or -v pool=PDOs_Treated"
  exit 1
fi

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
cd "$WD"
Rscript analysis/demultiplex/Auto_04_write_demultiplexed_counts.R --pool "$pool"

echo $(date +%T)
