#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=08:00:00
#PBS -N Auto_PDO_NBImg
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
OUT="${WD}/PDOs_outs/Auto_PDO_numbat"
SIF="${OUT}/Auto_numbat-rbase_latest.sif"
RLIB="${OUT}/Rlib"
NUMBAT_SRC="${OUT}/Auto_numbat_source"

mkdir -p "${OUT}/singularity_cache" "${OUT}/tmp" "${OUT}/logs" "$RLIB"
export SINGULARITY_CACHEDIR="${OUT}/singularity_cache"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export TMPDIR="${OUT}/tmp"
export APPTAINER_TMPDIR="${TMPDIR}"

if [[ ! -f "$SIF" ]]; then
  apptainer build --mksquashfs-args "-no-xattrs" "$SIF" docker://pkharchenkolab/numbat-rbase:latest
fi

if [[ ! -f "${NUMBAT_SRC}/DESCRIPTION" ]]; then
  mkdir -p "$NUMBAT_SRC"
  apptainer exec --cleanenv -B /rds:/rds "$SIF" bash -lc "cd /numbat && tar cf - ." | tar xf - -C "$NUMBAT_SRC"
fi

if ! apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" Rscript -e 'quit(status = ifelse(requireNamespace("numbat", quietly=TRUE), 0, 1))'; then
  apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" R CMD INSTALL -l "$RLIB" "$NUMBAT_SRC"
fi

apptainer exec --cleanenv --env R_LIBS_USER="$RLIB" -B /rds:/rds "$SIF" Rscript -e 'stopifnot(requireNamespace("numbat", quietly=TRUE)); library(numbat); cat(as.character(packageVersion("numbat")), "\n"); stopifnot(file.exists("/numbat/inst/bin/pileup_and_phase.R")); stopifnot(exists("run_numbat")); stopifnot(exists("ref_hca"))'

echo $(date +%T)
