#!/bin/bash
####################
# Analysis registry:
#   Status: active execution/support wrapper
#   Script: analysis/cnv/wes_subclone/Auto_02_submit_phylowgs.sh
#   Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
####################
# Submit PhyloWGS environment preparation and high-confidence sample runs.
####################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

mkdir -p "${OUT_ROOT}/logs"
cd "$PROJECT_DIR"

env_job="$(qsub -N Auto_PhyloEnv "${SCRIPT_DIR}/Auto_prepare_phylowgs_env.sh")"

for sample_name in "${HIGH_CONFIDENCE_SAMPLES[@]}"; do
  while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do
    sleep 180
  done
  qsub -W "depend=afterok:${env_job}" -v sample="$sample_name" -N "Phylo_${sample_name%%_vs_*}" "${SCRIPT_DIR}/Auto_run_phylowgs_sample.sh"
done
