#!/bin/bash
####################
# Analysis registry:
#   Status: legacy wrapper; retained for provenance, no current downstream use
#   Script: analysis/cnv/legacy_Auto_run_pdo_numbat_summary_plots.sh
#   Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=2:mem=32gb
#PBS -l walltime=04:00:00
#PBS -N Auto_PDO_NBSummary
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
cd "$WD"
Rscript analysis/cnv/legacy_Auto_PDO_numbat_concordance_summary_plots.R

echo $(date +%T)
