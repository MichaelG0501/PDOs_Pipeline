#!/bin/bash
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=24:00:00
#PBS -N Auto_WESFACETS
#PBS -koed
####################
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_run_wes_facets_sample.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - Sarek recalibrated tumour/normal CRAMs
#   - Common SNP VCF from Auto_prepare_facets_resources.sh
#   - Reference FASTA validated by Auto_validate_sarek_reference.sh
# Outputs:
#   - PDOs_outs/Auto_wes_subclone/intermediate/facets/<sample>/<sample>.snp_pileup.gz
#   - PDOs_outs/Auto_wes_subclone/tables/facets/Auto_<sample>_facets_*.tsv
# Downstream use: allele-specific FACETS/purity/ploidy input for conditional-shift WES visualization and HATCHet/THetA clone-CNA audits.
####################

set -euo pipefail

echo $(date +%T)

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>" >&2
  exit 1
fi

if ! wes_subclone_is_high_confidence "$sample"; then
  if wes_subclone_is_low_confidence "$sample" && [[ "${ALLOW_LOW_CONFIDENCE:-0}" != "1" ]]; then
    echo "ERROR: ${sample} is intentionally excluded from the default WES FACETS pipeline due to low PASS SNV support." >&2
    exit 1
  fi
  if ! wes_subclone_pair_for_sample "$sample" >/dev/null; then
    echo "ERROR: unknown sample: $sample" >&2
    exit 1
  fi
fi

if [[ -z "$WES_SUBCLONE_FACETS_SNP_VCF" || ! -f "$WES_SUBCLONE_FACETS_SNP_VCF" ]]; then
  echo "ERROR: WES_SUBCLONE_FACETS_SNP_VCF must point to an indexed common-SNP VCF for FACETS snp-pileup." >&2
  exit 1
fi
if [[ ! -f "${WES_SUBCLONE_FACETS_SNP_VCF}.tbi" && ! -f "${WES_SUBCLONE_FACETS_SNP_VCF}.csi" ]]; then
  echo "ERROR: FACETS SNP VCF is not indexed: $WES_SUBCLONE_FACETS_SNP_VCF" >&2
  exit 1
fi

module purge
module load tools/dev
module load tools/prod
module load BCFtools/1.22-GCC-14.2.0
module load SAMtools/1.22.1-GCC-14.2.0

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate "$WES_SUBCLONE_CONDA_ENV"

command -v "$SNP_PILEUP_BIN" >/dev/null || { echo "ERROR: snp-pileup not found. Set SNP_PILEUP_BIN or use an env containing FACETS snp-pileup." >&2; exit 1; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found after module load." >&2; exit 1; }
Rscript -e 'suppressPackageStartupMessages(library(facets)); suppressPackageStartupMessages(library(data.table))' >/dev/null

read -r tumor normal tumor_vcf_sample normal_vcf_sample tumor_cram normal_cram mutect_vcf cnvkit_cns < <(wes_subclone_paths_for_sample "$sample")
for path in "$tumor_cram" "$normal_cram"; do
  [[ -f "$path" ]] || { echo "ERROR: missing input: $path" >&2; exit 1; }
done

mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/intermediate/facets/${sample}" "${OUT_ROOT}/reports/${sample}" "${OUT_ROOT}/tmp"
export TMPDIR="${OUT_ROOT}/tmp"
export OUT_ROOT FACETS_CVAL FACETS_MIN_NHET EXCLUDE_SEX_CHROMS_DEFAULT
export REF_PATH="$WES_SUBCLONE_REF_FASTA"
export REF_CACHE="$WES_SUBCLONE_REF_CACHE"

"${SCRIPT_DIR}/Auto_validate_sarek_reference.sh" "$sample"

pileup="${OUT_ROOT}/intermediate/facets/${sample}/${sample}.snp_pileup.gz"
if [[ "${PDO_FORCE_REBUILD:-0}" == "1" || ! -s "$pileup" ]]; then
  echo "Running FACETS snp-pileup for ${sample}"
  "$SNP_PILEUP_BIN" ${FACETS_SNP_PILEUP_OPTS} \
    "$WES_SUBCLONE_FACETS_SNP_VCF" \
    "$pileup" \
    "$normal_cram" \
    "$tumor_cram"
else
  echo "Reusing existing snp-pileup: $pileup"
fi
if [[ "$(gzip -dc "$pileup" | wc -l)" -le 1 ]]; then
  echo "ERROR: snp-pileup wrote no SNP count rows for ${sample}. Check CRAM reference cache and SNP VCF compatibility." >&2
  exit 1
fi

Rscript "${SCRIPT_DIR}/Auto_run_facets_sample.R" "$sample" "$OUT_ROOT" "$FACETS_CVAL"

{
  printf "sample\t%s\n" "$sample"
  printf "tumor_cram\t%s\n" "$tumor_cram"
  printf "normal_cram\t%s\n" "$normal_cram"
  printf "facets_snp_vcf\t%s\n" "$WES_SUBCLONE_FACETS_SNP_VCF"
  printf "reference_fasta\t%s\n" "$WES_SUBCLONE_REF_FASTA"
  printf "facets_segments\t%s\n" "${OUT_ROOT}/tables/facets/Auto_${sample}_facets_segments.tsv"
  printf "facets_purity_ploidy\t%s\n" "${OUT_ROOT}/tables/facets/Auto_${sample}_facets_purity_ploidy.tsv"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${OUT_ROOT}/logs/Auto_${sample}_wes_facets_done.tsv"

echo $(date +%T)
