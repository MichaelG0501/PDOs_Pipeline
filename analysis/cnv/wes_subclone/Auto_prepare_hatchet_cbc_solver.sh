#!/bin/bash
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=02:00:00
#PBS -N Auto_HATCHetCBC
#PBS -koed
####################
# Analysis registry
# Status: active setup
# Script: analysis/cnv/wes_subclone/Auto_prepare_hatchet_cbc_solver.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - existing ephemeral HATCHet conda environment
#   - conda-forge package coin-or-cbc
# Outputs:
#   - CBC executable in ephemeral HATCHet conda environment
#   - live PDOs_outs/Auto_wes_clone_cna/logs/Auto_hatchet_cbc_solver.tsv
# Downstream use: HATCHet compute-cn with Pyomo/CBC instead of broken C++/Gurobi solver.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

LIVE_PROJECT="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
EPHEMERAL_PROJECT="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
LIVE_OUT="${LIVE_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
EPHEMERAL_OUT="${EPHEMERAL_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
HATCHET_ENV_DIR="${EPHEMERAL_OUT}/hatchet_env"
LOG_DIR="${LIVE_OUT}/logs"
STATUS_FILE="${LOG_DIR}/Auto_hatchet_cbc_solver.tsv"

mkdir -p "$LOG_DIR"

eval "$(~/miniforge3/bin/conda shell.bash hook)"
if [[ ! -x "${HATCHET_ENV_DIR}/bin/hatchet" ]]; then
  {
    printf "field\tvalue\n"
    printf "status\tblocked_missing_hatchet_env\n"
    printf "hatchet_env_dir\t%s\n" "$HATCHET_ENV_DIR"
    printf "finished\t%s\n" "$(date -Iseconds)"
  } > "$STATUS_FILE"
  cat "$STATUS_FILE"
  exit 1
fi

if [[ ! -x "${HATCHET_ENV_DIR}/bin/cbc" ]]; then
  conda install -y -p "$HATCHET_ENV_DIR" -c conda-forge coin-or-cbc
fi

source activate "$HATCHET_ENV_DIR"
{
  printf "field\tvalue\n"
  printf "status\tcomplete\n"
  printf "hatchet_env_dir\t%s\n" "$HATCHET_ENV_DIR"
  printf "cbc_path\t%s\n" "$(command -v cbc || true)"
  printf "cbc_version\t%s\n" "$(cbc -stop 2>&1 | head -1 || true)"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "$STATUS_FILE"
conda deactivate

cat "$STATUS_FILE"
echo $(date +%T)
