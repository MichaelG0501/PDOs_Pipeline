#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=04:00:00
#PBS -N Auto_ThetaN3Plot
#PBS -koed
####################
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_run_theta2_n3_clone_cna_compare.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - analysis/cnv/wes_subclone/Auto_theta2_n3_clone_cna_compare.R
#   - PDOs_outs/Auto_wes_clone_cna/tables/Auto_theta2_n3_clone_cna_manifest.csv
# Outputs:
#   - PDOs_outs/Auto_wes_clone_cna/figures/Auto_theta2_n3_clone_cna_compare_*.pdf/.png
#   - PDOs_outs/Auto_wes_clone_cna/tables/Auto_theta2_n3_clone_cna_*.csv
# Downstream use: terminal visual assessment of forced THetA2 n3 clone-CNA profiles.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd "$WD"
Rscript analysis/cnv/wes_subclone/Auto_theta2_n3_clone_cna_compare.R
echo $(date +%T)
