#!/bin/bash
#PBS -l select=1:ncpus=4:mem=24gb
#PBS -l walltime=04:00:00
#PBS -N Auto_ThetaEnv
#PBS -koed
####################
# Analysis registry
# Status: active setup
# Script: analysis/cnv/wes_subclone/Auto_prepare_theta2_clone_cna_env.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs: Bioconda/conda-forge package repositories.
# Outputs:
#   - ephemeral PDOs_outs/Auto_wes_clone_cna/theta2_env
#   - ephemeral PDOs_outs/Auto_wes_clone_cna/cnvkit_py310_env
#   - live PDOs_outs/Auto_wes_clone_cna/logs/Auto_theta2_clone_cna_env.tsv
# Downstream use: software setup for THetA2/CNVkit clone-specific WES CNA.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"

LIVE_PROJECT="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
EPHEMERAL_PROJECT="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
LIVE_OUT="${LIVE_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
EPHEMERAL_OUT="${EPHEMERAL_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
THETA_ENV_DIR="${EPHEMERAL_OUT}/theta2_env"
CNVKIT_ENV_DIR="${EPHEMERAL_OUT}/cnvkit_py310_env"
LOG_DIR="${LIVE_OUT}/logs"
LOG_FILE="${LOG_DIR}/Auto_theta2_clone_cna_env.tsv"

mkdir -p "$LOG_DIR" "$EPHEMERAL_OUT"

write_status() {
  local status="$1"
  local message="$2"
  {
    printf "field\tvalue\n"
    printf "status\t%s\n" "$status"
    printf "message\t%s\n" "$message"
    printf "theta_env_dir\t%s\n" "$THETA_ENV_DIR"
    printf "cnvkit_env_dir\t%s\n" "$CNVKIT_ENV_DIR"
    printf "cnvkit\t%s\n" "$(command -v cnvkit.py 2>/dev/null || true)"
    printf "RunTHetA.py\t%s\n" "$(command -v RunTHetA.py 2>/dev/null || true)"
    printf "python\t%s\n" "$(command -v python 2>/dev/null || true)"
    printf "finished\t%s\n" "$(date -Iseconds)"
  } > "$LOG_FILE"
}

if [[ ! -x "${THETA_ENV_DIR}/bin/RunTHetA.py" ]]; then
  rm -f "$LOG_FILE"
  mamba create -y -p "$THETA_ENV_DIR" \
    -c conda-forge -c bioconda --strict-channel-priority \
    "theta2"
fi

set +u
source activate "$THETA_ENV_DIR"
set -u
command -v RunTHetA.py >/dev/null
RunTHetA.py --help >/dev/null 2>&1 || true
conda deactivate

if [[ ! -x "${CNVKIT_ENV_DIR}/bin/cnvkit.py" ]]; then
  mamba create -y -p "$CNVKIT_ENV_DIR" \
    -c conda-forge -c bioconda --strict-channel-priority \
    "python=3.10" \
    "pandas<2" \
    "numpy<2" \
    "cnvkit"
fi

set +u
source activate "$CNVKIT_ENV_DIR"
set -u
command -v cnvkit.py >/dev/null
cnvkit.py version
conda deactivate

write_status "ready" "THetA2/CNVkit environment is available."
cat "$LOG_FILE"
echo $(date +%T)
