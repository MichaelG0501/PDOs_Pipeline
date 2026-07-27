#!/bin/bash
####################
# Create the local conda environment required by the WES subclone workflow.
#
# The environment is installed under PDOs_outs/Auto_wes_subclone/conda_env so
# no shared conda environment is modified.
####################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

env_dir="$WES_SUBCLONE_CONDA_ENV"
mkdir -p "$OUT_ROOT" "${OUT_ROOT}/logs"

eval "$(~/miniforge3/bin/conda shell.bash hook)"

if [[ ! -x "${env_dir}/bin/Rscript" || ! -x "${env_dir}/bin/pyclone-vi" ]]; then
  echo "Creating WES subclone environment: ${env_dir}"
  mamba create -y -p "$env_dir" \
    -c conda-forge -c bioconda --strict-channel-priority \
    "r-base=4.4" \
    "r-data.table" \
    "r-facets=0.6.2" \
    "pyclone-vi=0.2.0" \
    "bcftools=1.22" \
    "samtools=1.22" \
    "htslib"
else
  echo "Reusing existing WES subclone environment: ${env_dir}"
fi

source activate "$env_dir"

Rscript --vanilla -e 'suppressPackageStartupMessages(library(facets)); suppressPackageStartupMessages(library(data.table)); cat("facets/data.table OK\n")'
command -v pyclone-vi >/dev/null || { echo "ERROR: pyclone-vi is not available in ${env_dir}" >&2; exit 1; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools is not available in ${env_dir}" >&2; exit 1; }
command -v samtools >/dev/null || { echo "ERROR: samtools is not available in ${env_dir}" >&2; exit 1; }

if ! command -v snp-pileup >/dev/null; then
  facets_snp_pileup="$(Rscript --vanilla -e 'p <- system.file("extcode", "snp-pileup", package = "facets"); if (nzchar(p)) cat(p)' 2>/dev/null)"
  if [[ -n "$facets_snp_pileup" && -x "$facets_snp_pileup" ]]; then
    ln -sf "$facets_snp_pileup" "${env_dir}/bin/snp-pileup"
  fi
fi
if ! command -v snp-pileup >/dev/null; then
  facets_extcode="$(Rscript --vanilla -e 'cat(system.file("extcode", package = "facets"))')"
  if [[ -f "${facets_extcode}/snp-pileup.cpp" ]]; then
    echo "Compiling FACETS snp-pileup from installed package source"
    "${CXX:-g++}" -std=c++11 \
      -I"${CONDA_PREFIX}/include" \
      "${facets_extcode}/snp-pileup.cpp" \
      -L"${CONDA_PREFIX}/lib" \
      -lhts \
      -Wl,-rpath,"${CONDA_PREFIX}/lib" \
      -o "${env_dir}/bin/snp-pileup"
  fi
fi
command -v snp-pileup >/dev/null || { echo "ERROR: snp-pileup is not available after installing r-facets." >&2; exit 1; }

{
  printf "resource\tpath_or_version\n"
  printf "WES_SUBCLONE_CONDA_ENV\t%s\n" "$env_dir"
  printf "Rscript\t%s\n" "$(command -v Rscript)"
  printf "snp-pileup\t%s\n" "$(command -v snp-pileup)"
  printf "pyclone-vi\t%s\n" "$(command -v pyclone-vi)"
  printf "bcftools\t%s\n" "$(command -v bcftools)"
  printf "samtools\t%s\n" "$(command -v samtools)"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${OUT_ROOT}/logs/Auto_wes_subclone_env.tsv"

cat "${OUT_ROOT}/logs/Auto_wes_subclone_env.tsv"
