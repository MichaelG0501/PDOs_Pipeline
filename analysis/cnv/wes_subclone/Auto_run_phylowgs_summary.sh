#!/bin/bash
####################
# Analysis registry:
#   Status: active execution/support wrapper
#   Script: analysis/cnv/wes_subclone/Auto_run_phylowgs_summary.sh
#   Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -l walltime=01:00:00
#PBS -N Auto_PhyloSum
#PBS -koed
####################
# Summarise completed PhyloWGS result bundles into compact live tables.
####################

set -euo pipefail

echo $(date +%T)

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

cd "$PROJECT_DIR"
Rscript "${SCRIPT_DIR}/Auto_summarise_phylowgs_results.R" "$OUT_ROOT"

echo $(date +%T)
