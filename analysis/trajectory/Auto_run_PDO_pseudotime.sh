#!/bin/bash
#PBS -l select=1:ncpus=8:mem=96gb
#PBS -l walltime=24:00:00
#PBS -N Auto_PDO_pseudotime
#PBS -koed

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD
Rscript analysis/trajectory/Auto_PDO_pseudotime_samples.R
Rscript analysis/trajectory/Auto_PDO_pseudotime_linear_plot.R
Rscript analysis/trajectory/Auto_PDO_pseudotime_state_distance_matrix.R
echo $(date +%T)
