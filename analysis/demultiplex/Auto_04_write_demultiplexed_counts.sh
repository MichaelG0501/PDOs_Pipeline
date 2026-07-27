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

####################
# Scripts and durable manifests live in the live project; large matrices are
# read from the explicit ephemeral path by the R script.
WD="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
temporary_cluster_prefix="${temporary_cluster_prefix:-}"
counts_out_dir="${counts_out_dir:-}"
write_args=(--pool "$pool")
if [[ -n "$temporary_cluster_prefix" ]]; then
  write_args+=(--temporary_cluster_prefix "$temporary_cluster_prefix")
fi
if [[ -n "$counts_out_dir" ]]; then
  write_args+=(--counts_out_dir "$counts_out_dir")
fi
####################
cd "$WD"
Rscript analysis/demultiplex/Auto_04_write_demultiplexed_counts.R "${write_args[@]}"

echo $(date +%T)
