#!/bin/bash
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=02:00:00
#PBS -N Auto_SubMat
#PBS -koed
####################
# PBS wrapper for WES <-> scRNA subclone matching visualisation.
# Reads large Numbat bulk_clones files from ephemeral, writes
# figures and tables to live.
####################

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd $WD
Rscript analysis/cnv/wes_subclone/legacy_Auto_wes_scrna_subclone_match.R
echo $(date +%T)
