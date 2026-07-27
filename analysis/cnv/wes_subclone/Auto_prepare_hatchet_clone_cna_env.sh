#!/bin/bash
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=04:00:00
#PBS -N Auto_HATCHetEnv
#PBS -koed
####################
# Analysis registry
# Status: active setup
# Script: analysis/cnv/wes_subclone/Auto_prepare_hatchet_clone_cna_env.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - Bioconda/conda-forge packages hatchet, samtools, bcftools, mosdepth, htslib
# Outputs:
#   - ephemeral PDOs_outs/Auto_wes_clone_cna/hatchet_env
#   - live PDOs_outs/Auto_wes_clone_cna/logs/Auto_hatchet_clone_cna_env.tsv
# Downstream use: HATCHet/HATCHet2 clone-specific CNA inference from WES data.
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
STATUS_FILE="${LOG_DIR}/Auto_hatchet_clone_cna_env.tsv"

mkdir -p "$LOG_DIR" "$EPHEMERAL_OUT"

eval "$(~/miniforge3/bin/conda shell.bash hook)"

if [[ ! -x "${HATCHET_ENV_DIR}/bin/hatchet" ]]; then
  conda create -y -p "$HATCHET_ENV_DIR" \
    -c conda-forge -c bioconda \
    "hatchet=2.2.0" samtools bcftools mosdepth htslib tabix
fi

source activate "$HATCHET_ENV_DIR"
{
  printf "field\tvalue\n"
  printf "status\tcomplete\n"
  printf "hatchet_env_dir\t%s\n" "$HATCHET_ENV_DIR"
  printf "hatchet_version\t%s\n" "$(hatchet --version 2>&1 || true)"
  printf "hatchet_path\t%s\n" "$(command -v hatchet || true)"
  printf "samtools_version\t%s\n" "$(samtools --version | head -1 || true)"
  printf "bcftools_version\t%s\n" "$(bcftools --version | head -1 || true)"
  printf "mosdepth_version\t%s\n" "$(mosdepth --version 2>&1 | head -1 || true)"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "$STATUS_FILE"
conda deactivate

cat "$STATUS_FILE"
echo $(date +%T)
