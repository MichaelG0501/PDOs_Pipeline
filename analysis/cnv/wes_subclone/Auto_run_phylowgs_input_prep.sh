#!/bin/bash
####################
# Analysis registry:
#   Status: active execution/support wrapper
#   Script: analysis/cnv/wes_subclone/Auto_run_phylowgs_input_prep.sh
#   Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=02:00:00
#PBS -N Auto_PhyloIn
#PBS -koed
####################
# Per-sample PhyloWGS input-preparation PBS wrapper.
#
# This is intentionally independent of the Python 2 PhyloWGS environment so
# the FACETS/PyClone handoff tables are produced even if tool installation is
# delayed.
####################

set -euo pipefail

echo $(date +%T)

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>" >&2
  exit 1
fi
if ! wes_subclone_is_high_confidence "$sample"; then
  echo "ERROR: ${sample} is not in the high-confidence WES subclone allow-list." >&2
  exit 1
fi

EPHEMERAL_ROOT="${EPHEMERAL_OUT_ROOT:-/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone}"

module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

cd "$PROJECT_DIR"
Rscript "${SCRIPT_DIR}/Auto_make_phylowgs_inputs.R" "$sample" "$OUT_ROOT" "$EPHEMERAL_ROOT"

echo $(date +%T)
