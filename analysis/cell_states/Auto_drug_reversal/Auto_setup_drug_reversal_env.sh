#!/bin/bash
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=04:00:00
#PBS -N Auto_Drug_Env
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD
OUT_DIR=$WD/PDOs_outs/Auto_drug_reversal
ENV_PREFIX=$OUT_DIR/conda/Auto_drug_reversal
mkdir -p "$OUT_DIR/conda"
mkdir -p "$OUT_DIR/resources"
if [[ ! -d "$ENV_PREFIX" ]]; then
  conda env create --prefix "$ENV_PREFIX" --file analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_environment.yml
else
  conda env update --prefix "$ENV_PREFIX" --file analysis/cell_states/Auto_drug_reversal/Auto_drug_reversal_environment.yml --prune
fi
source activate "$ENV_PREFIX"
Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org"); remotes::install_github("immunogenomics/presto", upgrade="never"); remotes::install_github("lanagarmire/Asgard", upgrade="never"); remotes::install_github("SDTC-CPMed/scDrugPrio", upgrade="never")'
curl -L --fail --show-error --output "$OUT_DIR/resources/Auto_lit_ppi.rda" https://raw.githubusercontent.com/SDTC-CPMed/scDrugPrio/main/data/lit_ppi.rda
curl -L --fail --show-error --output "$OUT_DIR/resources/Auto_drug_bank_example_data.rda" https://raw.githubusercontent.com/SDTC-CPMed/scDrugPrio/main/data/drug_bank_example_data.rda
echo $(date +%T)
