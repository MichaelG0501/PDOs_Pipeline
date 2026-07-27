#!/bin/bash
#PBS -l select=1:ncpus=4:mem=12gb
#PBS -l walltime=04:00:00
#PBS -N Auto_PhyloEnv
#PBS -koed
####################
# Prepare an ephemeral PhyloWGS Python 2/C++ environment.
#
# PhyloWGS is an old Python 2 workflow. Keep the tool checkout and conda env
# under ephemeral storage; live outputs only record logs and final run products.
####################

set -euo pipefail

echo $(date +%T)

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

EPHEMERAL_ROOT="${EPHEMERAL_OUT_ROOT:-/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone}"
PHYLOWGS_ENV="${PHYLOWGS_ENV:-${EPHEMERAL_ROOT}/phylowgs_env}"
PHYLOWGS_DIR="${PHYLOWGS_DIR:-${EPHEMERAL_ROOT}/tools/phylowgs}"
LOG_DIR="${OUT_ROOT}/logs"
mkdir -p "$LOG_DIR" "$(dirname "$PHYLOWGS_DIR")"

module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"

if [[ ! -x "${PHYLOWGS_ENV}/bin/python2" ]]; then
  mamba create -y -p "$PHYLOWGS_ENV" -c conda-forge python=2.7 numpy scipy gsl ete2 pip
fi

conda activate "$PHYLOWGS_ENV"
export LD_LIBRARY_PATH="${PHYLOWGS_ENV}/lib:${LD_LIBRARY_PATH:-}"

if [[ ! -d "$PHYLOWGS_DIR/.git" ]]; then
  git clone https://github.com/morrislab/phylowgs.git "$PHYLOWGS_DIR"
else
  git -C "$PHYLOWGS_DIR" fetch --all --tags
fi

cd "$PHYLOWGS_DIR"
g++ -o mh.o -O3 mh.cpp util.cpp $(gsl-config --cflags --libs)
ldd mh.o
python2 - <<'PY'
import numpy
import scipy
import ete2
print("python2_ok")
PY

{
  printf "phylowgs_env\t%s\n" "$PHYLOWGS_ENV"
  printf "phylowgs_dir\t%s\n" "$PHYLOWGS_DIR"
  printf "mh_binary\t%s\n" "${PHYLOWGS_DIR}/mh.o"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${LOG_DIR}/Auto_phylowgs_env.tsv"

echo $(date +%T)
