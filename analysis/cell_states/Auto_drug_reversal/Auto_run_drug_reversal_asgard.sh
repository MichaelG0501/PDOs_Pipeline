#!/bin/bash
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -l walltime=24:00:00
#PBS -N Auto_ASGARD
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
if [[ -f PDOs_outs/Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.sh ]]; then
  source PDOs_outs/Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.sh
fi
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_asgard.R
echo $(date +%T)
