#!/bin/bash
####################
# Analysis registry:
#   Status: legacy wrapper; retained for provenance, no current downstream use
#   Script: analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_local_cmap.sh
#   Methodology: analysis/methodology/cell_states/Auto_drug_reversal_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=2:mem=64gb
#PBS -l walltime=04:00:00
#PBS -N Auto_Drug_CMap
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
ENV_PREFIX=/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Auto_drug_reversal/conda/Auto_drug_reversal
if [[ -d "$ENV_PREFIX" ]]; then
  source activate "$ENV_PREFIX"
else
  source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
fi
cd $WD
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_local_cmap.R
echo $(date +%T)
