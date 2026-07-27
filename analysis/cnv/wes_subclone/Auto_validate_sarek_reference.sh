#!/bin/bash
####################
# Validate that the candidate GRCh38 FASTA is compatible with the Sarek CRAMs.
#
# The original Sarek launch is parsed and recorded before any reference checks.
# FACETS wrappers call this script before snp-pileup. It exits non-zero if the
# FASTA does not match CRAM @SQ SN/LN/M5 fields or cannot decode a CRAM slice.
####################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

sample="${1:-${sample:-}}"
if [[ -z "$sample" ]]; then
  echo "ERROR: supply sample as first argument or environment variable sample" >&2
  exit 1
fi

if ! wes_subclone_pair_for_sample "$sample" >/dev/null; then
  echo "ERROR: unknown sample: $sample" >&2
  exit 1
fi

if [[ -z "$WES_SUBCLONE_REF_FASTA" ]]; then
  echo "ERROR: WES_SUBCLONE_REF_FASTA is not set. Provide the exact GATK.GRCh38-compatible FASTA used for CRAM decoding." >&2
  exit 1
fi
if [[ ! -f "$WES_SUBCLONE_REF_FASTA" ]]; then
  echo "ERROR: reference FASTA does not exist: $WES_SUBCLONE_REF_FASTA" >&2
  exit 1
fi
if [[ ! -f "${WES_SUBCLONE_REF_FASTA}.fai" ]]; then
  echo "ERROR: reference FASTA index is missing: ${WES_SUBCLONE_REF_FASTA}.fai" >&2
  echo "       Create/use an indexed reference in a writable reference location before running FACETS." >&2
  exit 1
fi
if [[ ! -f "$SAREK_LAUNCH" ]]; then
  echo "ERROR: original Sarek launch script is missing: $SAREK_LAUNCH" >&2
  exit 1
fi

read -r tumor normal tumor_vcf_sample normal_vcf_sample tumor_cram normal_cram mutect_vcf cnvkit_cns < <(wes_subclone_paths_for_sample "$sample")
for path in "$tumor_cram" "${tumor_cram}.crai" "$normal_cram" "${normal_cram}.crai"; do
  if [[ ! -f "$path" ]]; then
    echo "ERROR: missing CRAM/CRAI input: $path" >&2
    exit 1
  fi
done

val_dir="${OUT_ROOT}/intermediate/reference_validation/${sample}"
mkdir -p "$val_dir"

module purge
module load tools/prod
module load SAMtools/1.22.1-GCC-14.2.0

launch_flat="${val_dir}/Auto_${sample}_sarek_launch_command.txt"
tr '\n' ' ' < "$SAREK_LAUNCH" | sed 's/[[:space:]][[:space:]]*/ /g' > "$launch_flat"

extract_flag() {
  local flag="$1"
  awk -v flag="$flag" '{
    for (i = 1; i <= NF; i++) {
      if ($i == flag && i < NF) {
        print $(i + 1)
        found = 1
      } else if ($i ~ "^" flag "=") {
        sub("^" flag "=", "", $i)
        print $i
        found = 1
      }
    }
    if (!found) print "ABSENT"
  }' "$launch_flat"
}

{
  printf "parameter\tvalue\n"
  printf "sarek_launch\t%s\n" "$SAREK_LAUNCH"
  printf "nextflow_revision\t%s\n" "$(awk '{for(i=1;i<=NF;i++) if($i=="-r" && i<NF) print $(i+1)}' "$launch_flat")"
  printf "profile\t%s\n" "$(awk '{for(i=1;i<=NF;i++) if($i=="-profile" && i<NF) print $(i+1)}' "$launch_flat")"
  printf "genome\t%s\n" "$(extract_flag "--genome")"
  printf "fasta\t%s\n" "$(extract_flag "--fasta")"
  printf "igenomes_base\t%s\n" "$(extract_flag "--igenomes_base")"
  printf "igenomes_ignore\t%s\n" "$(extract_flag "--igenomes_ignore")"
  printf "wes_flag_present\t%s\n" "$(grep -Eq '(^|[[:space:]])--WES([[:space:]]|$)|(^|[[:space:]])--wes([[:space:]]|$)' "$launch_flat" && echo TRUE || echo FALSE)"
  printf "tools\t%s\n" "$(extract_flag "--tools")"
  printf "candidate_reference\t%s\n" "$WES_SUBCLONE_REF_FASTA"
} > "${val_dir}/Auto_${sample}_sarek_reference_parameters.tsv"

samtools dict "$WES_SUBCLONE_REF_FASTA" |
  awk 'BEGIN{OFS="\t"} /^@SQ/ {
    sn = ln = m5 = ""
    for (i = 1; i <= NF; i++) {
      if ($i ~ /^SN:/) { sn = substr($i, 4) }
      if ($i ~ /^LN:/) { ln = substr($i, 4) }
      if ($i ~ /^M5:/) { m5 = substr($i, 4) }
    }
    if (sn != "") print sn, ln, m5
  }' > "${val_dir}/Auto_${sample}_reference_sq.tsv"

validate_one_cram() {
  local label="$1"
  local cram="$2"
  local cram_sq="${val_dir}/Auto_${sample}_${label}_cram_sq.tsv"
  local mismatch="${val_dir}/Auto_${sample}_${label}_reference_mismatches.tsv"

  samtools view -H "$cram" |
    awk 'BEGIN{OFS="\t"} /^@SQ/ {
      sn = ln = m5 = ""
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^SN:/) { sn = substr($i, 4) }
        if ($i ~ /^LN:/) { ln = substr($i, 4) }
        if ($i ~ /^M5:/) { m5 = substr($i, 4) }
      }
      if (sn != "") print sn, ln, m5
    }' > "$cram_sq"

  awk 'BEGIN{FS=OFS="\t"}
    NR == FNR { ref_ln[$1] = $2; ref_m5[$1] = $3; next }
    {
      status = "OK"
      if (!($1 in ref_ln)) {
        status = "MISSING_CONTIG"
      } else if (ref_ln[$1] != $2) {
        status = "LENGTH_MISMATCH"
      } else if ($3 != "" && ref_m5[$1] != "" && ref_m5[$1] != $3) {
        status = "MD5_MISMATCH"
      }
      if (status != "OK") print $1, status, "cram_LN=" $2, "ref_LN=" ref_ln[$1], "cram_M5=" $3, "ref_M5=" ref_m5[$1]
    }' "${val_dir}/Auto_${sample}_reference_sq.tsv" "$cram_sq" > "$mismatch"

  if [[ -s "$mismatch" ]]; then
    echo "ERROR: reference FASTA does not match $label CRAM @SQ fields. See: $mismatch" >&2
    exit 1
  fi

  samtools view -T "$WES_SUBCLONE_REF_FASTA" "$cram" chr1:962000-962020 >/dev/null
}

validate_one_cram "tumor" "$tumor_cram"
validate_one_cram "normal" "$normal_cram"

{
  printf "sample\t%s\n" "$sample"
  printf "tumor_cram\t%s\n" "$tumor_cram"
  printf "normal_cram\t%s\n" "$normal_cram"
  printf "reference_fasta\t%s\n" "$WES_SUBCLONE_REF_FASTA"
  printf "status\tvalidated\n"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${val_dir}/Auto_${sample}_reference_validation_done.tsv"

echo "Reference validation passed for ${sample}"
