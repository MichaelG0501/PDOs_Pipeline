#!/bin/bash
#PBS -l select=1:ncpus=4:mem=48gb
#PBS -l walltime=04:00:00
#PBS -N Auto_PhyloNumbat
#PBS -koed
####################
# Run PhyloWGS inherited clone-CNA versus native Numbat CNA visualisation.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd "$WD"

Rscript analysis/cnv/wes_subclone/Auto_plot_phylowgs_numbat_compare.R

echo $(date +%T)
