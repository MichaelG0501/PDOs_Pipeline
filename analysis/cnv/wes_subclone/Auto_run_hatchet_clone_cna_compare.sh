#!/bin/bash
#PBS -l select=1:ncpus=4:mem=24gb
#PBS -l walltime=04:00:00
#PBS -N Auto_HATCHet_cmp
#PBS -koed

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd $WD
Rscript analysis/cnv/wes_subclone/Auto_hatchet_clone_cna_compare.R
echo $(date +%T)
