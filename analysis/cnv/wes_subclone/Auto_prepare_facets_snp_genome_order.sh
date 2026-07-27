#!/bin/bash
#PBS -l select=1:ncpus=4:mem=24gb
#PBS -l walltime=04:00:00
#PBS -N Auto_FACETSSort
#PBS -koed
####################
# Analysis registry
# Status: active utility
# Script: analysis/cnv/wes_subclone/Auto_prepare_facets_snp_genome_order.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - PDOs_outs/Auto_wes_subclone/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.vcf.gz
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.genome_order.vcf.gz
# Downstream use: corrected FACETS snp-pileup input sorted in reference
#   contig order rather than lexicographic chromosome order.
####################

set -euo pipefail

echo $(date +%T)

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

module purge
module load tools/prod
module load BCFtools/1.22-GCC-14.2.0
module load SAMtools/1.22.1-GCC-14.2.0

EPHEMERAL_ROOT="${EPHEMERAL_OUT_ROOT:-/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone}"
in_vcf="${OUT_ROOT}/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.vcf.gz"
if [[ ! -s "$in_vcf" ]]; then
  in_vcf="${EPHEMERAL_ROOT}/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.vcf.gz"
fi
out_vcf="${OUT_ROOT}/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.genome_order.vcf.gz"
tmp_dir="${OUT_ROOT}/tmp/Auto_facets_snp_genome_order"
log_dir="${OUT_ROOT}/logs"
mkdir -p "$tmp_dir" "$log_dir" "$(dirname "$out_vcf")"

if [[ ! -s "$in_vcf" ]]; then
  echo "ERROR: missing input VCF: $in_vcf" >&2
  exit 1
fi

if [[ "${PDO_FORCE_REBUILD:-0}" == "1" || ! -s "$out_vcf" || ! -s "${out_vcf}.tbi" ]]; then
  bcftools sort -T "$tmp_dir" -Oz -o "${out_vcf}.tmp" "$in_vcf"
  mv "${out_vcf}.tmp" "$out_vcf"
  tabix -f -p vcf "$out_vcf"
fi

{
  printf "input_vcf\t%s\n" "$in_vcf"
  printf "output_vcf\t%s\n" "$out_vcf"
  printf "output_index\t%s\n" "${out_vcf}.tbi"
  printf "first_records\t"
  bcftools view -H "$out_vcf" | cut -f1 | head -n 12 | paste -sd "," -
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${log_dir}/Auto_facets_snp_genome_order.tsv"

echo $(date +%T)
