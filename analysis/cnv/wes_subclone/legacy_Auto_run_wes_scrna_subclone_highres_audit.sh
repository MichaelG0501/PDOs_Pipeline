#!/bin/bash
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
