#!/bin/bash
#PBS -l select=1:ncpus=8:mem=256gb
#PBS -l walltime=04:00:00
#PBS -N update_optimal_mp
#PBS -o /rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/temp/
#PBS -e /rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/temp/
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD
Rscript analysis/metaprograms/Auto_ucell_vlnplot.R
echo $(date +%T)
