#!/bin/bash
#PBS -l select=1:ncpus=16:mem=96gb
#PBS -l walltime=36:00:00
#PBS -N Auto_WESSub
#PBS -koed
####################
# Per-sample FACETS -> PyClone-VI PBS wrapper.
#
# Default execution is restricted to PDO_1090_vs_NT_1090 and
# PDO_1181_vs_NT_1181. Set ALLOW_LOW_CONFIDENCE=1 only for deliberate manual
# experiments on lower-confidence samples.
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
    echo "ERROR: ${sample} is intentionally excluded from the default WES subclone pipeline due to low PASS SNV support." >&2
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
command -v "$PYCLONE_VI_BIN" >/dev/null || { echo "ERROR: pyclone-vi not found. Set PYCLONE_VI_BIN or use an env containing PyClone-VI." >&2; exit 1; }
command -v bcftools >/dev/null || { echo "ERROR: bcftools not found after module load." >&2; exit 1; }
Rscript -e 'suppressPackageStartupMessages(library(facets)); suppressPackageStartupMessages(library(data.table))' >/dev/null

read -r tumor normal tumor_vcf_sample normal_vcf_sample tumor_cram normal_cram mutect_vcf cnvkit_cns < <(wes_subclone_paths_for_sample "$sample")
for path in "$tumor_cram" "$normal_cram" "$mutect_vcf"; do
  [[ -f "$path" ]] || { echo "ERROR: missing input: $path" >&2; exit 1; }
done

mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/intermediate/facets/${sample}" "${OUT_ROOT}/reports/${sample}" "${OUT_ROOT}/tmp"
export TMPDIR="${OUT_ROOT}/tmp"
export OUT_ROOT FACETS_CVAL FACETS_MIN_NHET EXCLUDE_SEX_CHROMS_DEFAULT PYCLONE_MIN_TUMOUR_DEPTH PYCLONE_MIN_ALT_COUNT
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
MUTECT_VCF="$mutect_vcf" Rscript "${SCRIPT_DIR}/legacy_Auto_make_pyclone_input.R" "$sample" "$OUT_ROOT" "$mutect_vcf"

pyclone_dir="${OUT_ROOT}/intermediate/pyclone/${sample}"
mkdir -p "$pyclone_dir" "${OUT_ROOT}/tables/pyclone" "${OUT_ROOT}/reports/${sample}"
pyclone_input="${OUT_ROOT}/tables/pyclone/Auto_${sample}_pyclone_vi_input.tsv"
pyclone_h5="${pyclone_dir}/Auto_${sample}_pyclone_vi.h5"
pyclone_results="${OUT_ROOT}/tables/pyclone/Auto_${sample}_pyclone_vi_results.tsv"

if [[ "${PDO_FORCE_REBUILD:-0}" == "1" || ! -s "$pyclone_results" ]]; then
  "$PYCLONE_VI_BIN" fit \
    -i "$pyclone_input" \
    -o "$pyclone_h5" \
    -c "$PYCLONE_NUM_CLUSTERS" \
    -d "$PYCLONE_DENSITY" \
    -r "$PYCLONE_RESTARTS"
  "$PYCLONE_VI_BIN" write-results-file \
    -i "$pyclone_h5" \
    -o "$pyclone_results"
else
  echo "Reusing existing PyClone-VI results: $pyclone_results"
fi

{
  printf "sample\t%s\n" "$sample"
  printf "tumor_cram\t%s\n" "$tumor_cram"
  printf "normal_cram\t%s\n" "$normal_cram"
  printf "mutect_vcf\t%s\n" "$mutect_vcf"
  printf "facets_snp_vcf\t%s\n" "$WES_SUBCLONE_FACETS_SNP_VCF"
  printf "reference_fasta\t%s\n" "$WES_SUBCLONE_REF_FASTA"
  printf "pyclone_input\t%s\n" "$pyclone_input"
  printf "pyclone_results\t%s\n" "$pyclone_results"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${OUT_ROOT}/logs/Auto_${sample}_wes_subclone_done.tsv"

echo $(date +%T)
