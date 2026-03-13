#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=04:00:00
#PBS -N enrichment_annotation
#PBS -o /rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/temp/
#PBS -e /rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/temp/
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD
Rscript analysis/enrichment/Auto_enrichment_annotation.R
echo $(date +%T)
