#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=12:00:00
#PBS -N Auto_PDO_NBSubMP
#PBS -koed

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD
Rscript analysis/cnv/Auto_PDO_numbat_subclone_mp_heatmap.R
echo $(date +%T)
