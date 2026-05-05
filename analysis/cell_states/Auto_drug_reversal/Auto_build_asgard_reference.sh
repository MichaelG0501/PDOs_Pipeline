#!/bin/bash
#PBS -l select=1:ncpus=4:mem=160gb
#PBS -l walltime=36:00:00
#PBS -N Auto_ASGARD_Ref
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
ENV_PREFIX=$WD/PDOs_outs/Auto_drug_reversal/conda/Auto_drug_reversal
if [[ -d "$ENV_PREFIX" ]]; then
  source activate "$ENV_PREFIX"
else
  source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
fi
cd $WD
Rscript analysis/cell_states/Auto_drug_reversal/Auto_prepare_asgard_reference.R
echo $(date +%T)
