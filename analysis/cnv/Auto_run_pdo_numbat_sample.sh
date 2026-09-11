#!/bin/bash
####################
# Analysis registry:
#   Status: active execution/support wrapper
#   Script: analysis/cnv/Auto_run_pdo_numbat_sample.sh
#   Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=12:mem=128gb
#PBS -l walltime=48:00:00
#PBS -N Auto_PDO_NBRun
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>"
  exit 1
fi

WD="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
OUT="${WD}/PDOs_outs/Auto_PDO_numbat"
MANIFEST="${OUT}/Auto_PDO_numbat_manifest.csv"
CACHE_OUT="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_PDO_numbat"
SIF="${CACHE_OUT}/Auto_numbat-rbase_latest.sif"
NCORES="${NCORES:-12}"
RLIB="${CACHE_OUT}/Rlib"

mkdir -p "${OUT}/logs" "${CACHE_OUT}/singularity_cache" "${CACHE_OUT}/tmp"
export SINGULARITY_CACHEDIR="${CACHE_OUT}/singularity_cache"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export TMPDIR="${CACHE_OUT}/tmp"

if [[ ! -f "$SIF" ]]; then
  echo "ERROR: missing Numbat container: $SIF"
  exit 1
fi

cd "$WD"
apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" \
  Rscript analysis/cnv/Auto_PDO_numbat_run_sample.R "$sample" "$MANIFEST" "$NCORES"

echo $(date +%T)
