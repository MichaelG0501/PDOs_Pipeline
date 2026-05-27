#!/bin/bash
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

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
OUT="${WD}/PDOs_outs/Auto_PDO_numbat"
MANIFEST="${OUT}/Auto_PDO_numbat_manifest.csv"
SIF="${OUT}/Auto_numbat-rbase_latest.sif"
NCORES="${NCORES:-12}"
RLIB="${OUT}/Rlib"

mkdir -p "${OUT}/logs" "${OUT}/singularity_cache" "${OUT}/tmp"
export SINGULARITY_CACHEDIR="${OUT}/singularity_cache"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export TMPDIR="${OUT}/tmp"

if [[ ! -f "$SIF" ]]; then
  echo "ERROR: missing Numbat container: $SIF"
  exit 1
fi

cd "$WD"
apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" \
  Rscript analysis/cnv/Auto_PDO_numbat_run_sample.R "$sample" "$MANIFEST" "$NCORES"

echo $(date +%T)
