#!/bin/bash
####################
# Analysis registry:
#   Status: active execution/support wrapper
#   Script: analysis/cnv/wes_subclone/Auto_run_wes_absolute_cna_compare.sh
#   Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=4:mem=48gb
#PBS -l walltime=08:00:00
#PBS -N Auto_WESAbsCNA
#PBS -koed
####################
# Run high-resolution ploidy-adjusted WES/scRNA CNA comparison.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd "$WD"

Rscript analysis/cnv/wes_subclone/Auto_wes_absolute_highres_subclone_compare.R

echo $(date +%T)
