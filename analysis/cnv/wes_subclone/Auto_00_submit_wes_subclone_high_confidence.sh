#!/bin/bash
####################
# Submit only the high-confidence PDO WES subclone samples.
#
# Required environment before running:
#   WES_SUBCLONE_REF_FASTA=/path/to/Homo_sapiens_assembly38.fasta
#   WES_SUBCLONE_FACETS_SNP_VCF=/path/to/common_snps.vcf.gz
####################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

if [[ -z "$WES_SUBCLONE_REF_FASTA" || ! -f "$WES_SUBCLONE_REF_FASTA" ]]; then
  echo "ERROR: export WES_SUBCLONE_REF_FASTA before submitting." >&2
  exit 1
fi
if [[ -z "$WES_SUBCLONE_FACETS_SNP_VCF" || ! -f "$WES_SUBCLONE_FACETS_SNP_VCF" ]]; then
  echo "ERROR: export WES_SUBCLONE_FACETS_SNP_VCF before submitting." >&2
  exit 1
fi

mkdir -p "${OUT_ROOT}/logs"
cd "$PROJECT_DIR"

for sample_name in "${HIGH_CONFIDENCE_SAMPLES[@]}"; do
  while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do
    sleep 180
  done
  qsub -V -v sample="$sample_name" -N "WESFACETS_${sample_name%%_vs_*}" "${SCRIPT_DIR}/Auto_run_wes_facets_sample.sh"
done
