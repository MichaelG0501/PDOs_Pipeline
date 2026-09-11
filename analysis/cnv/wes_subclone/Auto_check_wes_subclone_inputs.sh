#!/bin/bash
####################
# Analysis registry:
#   Status: active execution/support wrapper
#   Script: analysis/cnv/wes_subclone/Auto_check_wes_subclone_inputs.sh
#   Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
#   Map: analysis/ANALYSIS_MAP.md
#   Description: Orchestrates the command, environment, resources, and
#     dependencies documented below; it does not define new analytical logic.
####################
####################
# Read-only input inventory for the WES FACETS -> PyClone-VI workflow.
####################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

mkdir -p "${OUT_ROOT}/logs"
out="${OUT_ROOT}/logs/Auto_wes_subclone_input_inventory.tsv"
printf "sample\tdefault_status\ttumor\tnormal\tmutect_vcf\tmutect_tbi\ttumor_cram\ttumor_crai\tnormal_cram\tnormal_crai\tcnvkit_somatic_cns\n" > "$out"

for sample_name in "${HIGH_CONFIDENCE_SAMPLES[@]}" "${LOW_CONFIDENCE_SAMPLES[@]}"; do
  read -r tumor normal tumor_vcf_sample normal_vcf_sample tumor_cram normal_cram mutect_vcf cnvkit_cns < <(wes_subclone_paths_for_sample "$sample_name")
  status="excluded_low_confidence"
  wes_subclone_is_high_confidence "$sample_name" && status="included_high_confidence"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample_name" \
    "$status" \
    "$tumor" \
    "$normal" \
    "$([[ -f "$mutect_vcf" ]] && echo present || echo missing)" \
    "$([[ -f "${mutect_vcf}.tbi" || -f "${mutect_vcf}.csi" ]] && echo present || echo missing)" \
    "$([[ -f "$tumor_cram" ]] && echo present || echo missing)" \
    "$([[ -f "${tumor_cram}.crai" ]] && echo present || echo missing)" \
    "$([[ -f "$normal_cram" ]] && echo present || echo missing)" \
    "$([[ -f "${normal_cram}.crai" ]] && echo present || echo missing)" \
    "$([[ -f "$cnvkit_cns" ]] && echo present || echo missing)" >> "$out"
done

cat "$out"
