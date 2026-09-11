#!/bin/bash
####################
# Analysis registry:
#   Status: legacy wrapper; retained for provenance, no current downstream use
#   Script: analysis/cnv/legacy_Auto_run_pdo_numbat_subclone_mp.sh
#   Methodology: analysis/methodology/cnv/cnv_workflows_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=12:00:00
#PBS -N Auto_PDO_NBSubMP
#PBS -koed

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd $WD
Rscript analysis/cnv/legacy_Auto_PDO_numbat_subclone_mp_heatmap.R
echo $(date +%T)
