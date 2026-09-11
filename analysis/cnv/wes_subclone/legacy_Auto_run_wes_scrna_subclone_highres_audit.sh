#!/bin/bash
####################
# Analysis registry:
#   Status: legacy wrapper; retained for provenance, no current downstream use
#   Script: analysis/cnv/wes_subclone/legacy_Auto_run_wes_scrna_subclone_highres_audit.sh
#   Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=4:mem=48gb
#PBS -l walltime=08:00:00
#PBS -N Auto_WESscRNA_HR
#PBS -koed
####################
# Run high-resolution WES/scRNA CNA concordance audit.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd "$WD"

Rscript analysis/cnv/wes_subclone/legacy_Auto_wes_scrna_subclone_highres_audit.R

echo $(date +%T)
