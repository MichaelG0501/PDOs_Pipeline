#!/bin/bash
####################
# Shared configuration for the WES FACETS -> PyClone-VI workflow.
#
# Source this file from the submission and per-sample wrappers. Paths are kept
# here so the biological restrictions and sample allow-list stay consistent.
####################

PROJECT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
SAREK_LAUNCH="/rds/general/project/spatialtranscriptomics/live/WES_PDO/sarek/sarek.sh"
SAREK_ROOT="/rds/general/project/spatialtranscriptomics/live/sarek_mutect"
OUT_ROOT="${PROJECT_DIR}/PDOs_outs/Auto_wes_subclone"

HIGH_CONFIDENCE_SAMPLES=(
  "PDO_1090_vs_NT_1090"
  "PDO_1181_vs_NT_1181"
)

LOW_CONFIDENCE_SAMPLES=(
  "PDO_1070_vs_NT_1070"
  "PDO_1072_vs_NT_1072"
  "PDO_1121_vs_NT_1121"
  "PDO_1141_vs_NT_1141"
)

EXCLUDE_SEX_CHROMS_DEFAULT="${EXCLUDE_SEX_CHROMS_DEFAULT:-1}"

FACETS_CVAL="${FACETS_CVAL:-150}"
FACETS_MIN_NHET="${FACETS_MIN_NHET:-15}"
FACETS_SNP_PILEUP_OPTS="${FACETS_SNP_PILEUP_OPTS:--g -q15 -Q20 -P100 -r25,0}"

PYCLONE_MIN_TUMOUR_DEPTH="${PYCLONE_MIN_TUMOUR_DEPTH:-20}"
PYCLONE_MIN_ALT_COUNT="${PYCLONE_MIN_ALT_COUNT:-3}"
PYCLONE_NUM_CLUSTERS="${PYCLONE_NUM_CLUSTERS:-40}"
PYCLONE_RESTARTS="${PYCLONE_RESTARTS:-10}"
PYCLONE_DENSITY="${PYCLONE_DENSITY:-beta-binomial}"

WES_SUBCLONE_CONDA_ENV="${WES_SUBCLONE_CONDA_ENV:-${OUT_ROOT}/conda_env}"
WES_SUBCLONE_REF_FASTA="${WES_SUBCLONE_REF_FASTA:-${OUT_ROOT}/resources/reference/Homo_sapiens_assembly38.fasta}"
WES_SUBCLONE_REF_CACHE="${WES_SUBCLONE_REF_CACHE:-${OUT_ROOT}/resources/reference/ref_cache/%2s/%2s/%s}"
WES_SUBCLONE_FACETS_SNP_VCF="${WES_SUBCLONE_FACETS_SNP_VCF:-${OUT_ROOT}/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.vcf.gz}"
SNP_PILEUP_BIN="${SNP_PILEUP_BIN:-snp-pileup}"
PYCLONE_VI_BIN="${PYCLONE_VI_BIN:-pyclone-vi}"

wes_subclone_is_high_confidence() {
  local query="$1"
  local sample
  for sample in "${HIGH_CONFIDENCE_SAMPLES[@]}"; do
    [[ "$query" == "$sample" ]] && return 0
  done
  return 1
}

wes_subclone_is_low_confidence() {
  local query="$1"
  local sample
  for sample in "${LOW_CONFIDENCE_SAMPLES[@]}"; do
    [[ "$query" == "$sample" ]] && return 0
  done
  return 1
}

wes_subclone_pair_for_sample() {
  local sample="$1"
  case "$sample" in
    PDO_1090_vs_NT_1090) echo "PDO_1090 NT_1090 SUR1090_PDO_1090 SUR1090_NT_1090" ;;
    PDO_1181_vs_NT_1181) echo "PDO_1181 NT_1181 SUR1181_PDO_1181 SUR1181_NT_1181" ;;
    PDO_1070_vs_NT_1070) echo "PDO_1070 NT_1070 SUR1070_PDO_1070 SUR1070_NT_1070" ;;
    PDO_1072_vs_NT_1072) echo "PDO_1072 NT_1072 SUR1072_PDO_1072 SUR1072_NT_1072" ;;
    PDO_1121_vs_NT_1121) echo "PDO_1121 NT_1121 SUR1121_PDO_1121 SUR1121_NT_1121" ;;
    PDO_1141_vs_NT_1141) echo "PDO_1141 NT_1141 SUR1141_PDO_1141 SUR1141_NT_1141" ;;
    *) return 1 ;;
  esac
}

wes_subclone_paths_for_sample() {
  local sample="$1"
  local tumor normal tumor_vcf_sample normal_vcf_sample
  read -r tumor normal tumor_vcf_sample normal_vcf_sample < <(wes_subclone_pair_for_sample "$sample")
  local mutect_dir="${SAREK_ROOT}/variant_calling/mutect2/${sample}"
  local cnvkit_dir="${SAREK_ROOT}/variant_calling/cnvkit/${sample}"
  local tumor_cram="${SAREK_ROOT}/preprocessing/recalibrated/${tumor}/${tumor}.recal.cram"
  local normal_cram="${SAREK_ROOT}/preprocessing/recalibrated/${normal}/${normal}.recal.cram"
  local mutect_vcf="${mutect_dir}/${sample}.mutect2.filtered.vcf.gz"
  local cnvkit_cns="${cnvkit_dir}/${tumor}.somatic.call.cns"
  echo "$tumor $normal $tumor_vcf_sample $normal_vcf_sample $tumor_cram $normal_cram $mutect_vcf $cnvkit_cns"
}
