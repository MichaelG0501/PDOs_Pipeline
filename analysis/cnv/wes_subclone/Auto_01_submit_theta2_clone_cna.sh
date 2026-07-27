#!/bin/bash
####################
# Analysis registry
# Status: active submitter
# Script: analysis/cnv/wes_subclone/Auto_01_submit_theta2_clone_cna.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs: high-confidence WES subclone sample list from Auto_wes_subclone_config.sh.
# Outputs: PBS jobs for THetA2 clone-specific WES CNA.
# Downstream use: submission helper only.
####################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

cd "$PROJECT_DIR"

for sample_name in "${HIGH_CONFIDENCE_SAMPLES[@]}"; do
  while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do
    sleep 180
  done
  qsub -v sample="$sample_name" -N "Theta_${sample_name%%_vs_*}" "${SCRIPT_DIR}/Auto_run_theta2_clone_cna_sample.sh"
done
