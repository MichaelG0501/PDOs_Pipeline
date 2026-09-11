#!/bin/bash
####################
# Analysis registry:
#   Status: legacy wrapper; retained for provenance, no current downstream use
#   Script: analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_inputs.sh
#   Methodology: analysis/methodology/cell_states/Auto_drug_reversal_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=4:mem=96gb
#PBS -l walltime=08:00:00
#PBS -N Auto_Drug_Input
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
cd $WD
export AUTO_DRUG_DEG_MODE=pseudobulk
Rscript analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_inputs.R
echo $(date +%T)
