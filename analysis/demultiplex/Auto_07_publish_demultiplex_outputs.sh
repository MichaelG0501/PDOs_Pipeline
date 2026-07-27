#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=12:00:00
#PBS -N Auto_PDO_PublishDemux
#PBS -koed

####################
# Analysis registry
# Status: active durable-output publisher for completed multiplexed PDO pools.
# Script: analysis/demultiplex/Auto_07_publish_demultiplex_outputs.sh
# Methodology: analysis/methodology/demultiplex/demultiplex_methodology.md
# Map: analysis/ANALYSIS_MAP.md
# Inputs: completed Cell Ranger pool directory and exported demultiplex CSVs.
# Outputs: non-overwriting live copies in ITH_sc/PDOs/Cellranger_outs/<pool>/
#          and ITH_sc/PDOs/00_counts_matrix_all/.
# Downstream: both durable output locations are downstream analysis inputs.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

pool="${pool:-}"
cellranger_source="${cellranger_source:-}"
cellranger_destination_root="${cellranger_destination_root:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Cellranger_outs}"
counts_source_dir="${counts_source_dir:-}"
counts_destination_dir="${counts_destination_dir:-/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/00_counts_matrix_all}"

if [[ -z "$pool" || -z "$cellranger_source" || -z "$counts_source_dir" ]]; then
  echo "ERROR: submit with pool, cellranger_source, and counts_source_dir variables"
  exit 1
fi
if [[ ! -d "$cellranger_source" || ! -d "$counts_source_dir" ]]; then
  echo "ERROR: a required source directory is missing"
  exit 1
fi

####################
# Never overwrite or move existing data. The pool-specific destination matches
# the established historical Cellranger_outs layout.
cellranger_destination="${cellranger_destination_root}/${pool}"
if [[ -e "$cellranger_destination" ]]; then
  echo "ERROR: refusing to overwrite existing Cell Ranger destination: $cellranger_destination"
  exit 1
fi
mkdir -p "$cellranger_destination_root" "$counts_destination_dir"

mapfile -t count_files < <(find "$counts_source_dir" -maxdepth 1 -type f -name "TEMP_${pool}_SouporcellCluster*_PDO.csv" -printf '%f\n' | sort)
if [[ ${#count_files[@]} -eq 0 ]]; then
  echo "ERROR: no temporary-cluster count CSVs found in $counts_source_dir"
  exit 1
fi
for count_file in "${count_files[@]}"; do
  if [[ -e "${counts_destination_dir}/${count_file}" ]]; then
    echo "ERROR: refusing to overwrite existing count CSV: ${counts_destination_dir}/${count_file}"
    exit 1
  fi
done
####################

echo "Publishing Cell Ranger: $cellranger_source -> $cellranger_destination"
cp -a "$cellranger_source" "$cellranger_destination"
for count_file in "${count_files[@]}"; do
  echo "Publishing CSV: ${counts_source_dir}/${count_file} -> ${counts_destination_dir}/${count_file}"
  cp -a "${counts_source_dir}/${count_file}" "${counts_destination_dir}/${count_file}"
done
echo $(date +%T)
