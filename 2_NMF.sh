#!/bin/bash
#PBS -l select=1:ncpus=6:mem=40gb
#PBS -l walltime=4:00:00
#PBS -N nmf

echo $(date +%T)

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline

cd $WD

Rscript NMF.R $sample

echo $(date +%T)
