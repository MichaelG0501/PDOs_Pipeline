#!/bin/bash
####################
# Download and prepare reference resources required by the WES FACETS workflow.
#
# Outputs are written under PDOs_outs/Auto_wes_subclone/resources/. This script
# does not write into shared reference locations. It prepares:
#   1. Broad/GATK GRCh38 Homo_sapiens_assembly38 FASTA, FAI, and dict.
#   2. A bgzipped/tabix-indexed UCSC hg38 common biallelic SNP VCF for FACETS.
####################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

module purge
module load tools/prod
module load BCFtools/1.22-GCC-14.2.0
module load SAMtools/1.22.1-GCC-14.2.0

resource_root="${OUT_ROOT}/resources"
ref_dir="${resource_root}/reference"
snp_dir="${resource_root}/facets_snps"
log_dir="${OUT_ROOT}/logs"
mkdir -p "$ref_dir" "$snp_dir" "$log_dir"

ref_fasta="${ref_dir}/Homo_sapiens_assembly38.fasta"
ref_fai="${ref_fasta}.fai"
ref_dict="${ref_dir}/Homo_sapiens_assembly38.dict"
ref_gz="${ref_fasta}.gz"
ref_cache_dir="${ref_dir}/ref_cache"

gatk_ref_url="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz"
gatk_fai_url="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai"
gatk_dict_url="ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict"

download_if_missing() {
  local url="$1"
  local out="$2"
  if [[ -s "$out" ]]; then
    echo "Reusing existing: $out"
    return 0
  fi
  echo "Downloading: $url"
  curl -L --retry 3 --retry-delay 20 --fail -o "${out}.tmp" "$url"
  mv "${out}.tmp" "$out"
}

if [[ ! -s "$ref_fasta" ]]; then
  download_if_missing "$gatk_ref_url" "$ref_gz"
  echo "Decompressing reference FASTA"
  gzip -dc "$ref_gz" > "${ref_fasta}.tmp"
  mv "${ref_fasta}.tmp" "$ref_fasta"
fi
download_if_missing "$gatk_fai_url" "$ref_fai"
download_if_missing "$gatk_dict_url" "$ref_dict"
if [[ ! -d "$ref_cache_dir" || "$(find "$ref_cache_dir" -type f | wc -l)" -lt 24 ]]; then
  echo "Populating htslib REF_CACHE for CRAM decoding"
  mkdir -p "$ref_cache_dir"
  seq_cache_populate.pl -root "$ref_cache_dir" -subdirs 2 "$ref_fasta" > "${log_dir}/Auto_wes_subclone_ref_cache_populate.log"
fi

snp_txt="${snp_dir}/Auto_ucsc_hg38_snp151Common.txt.gz"
snp_sql="${snp_dir}/Auto_ucsc_hg38_snp151Common.sql"
snp_vcf="${snp_dir}/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.vcf.gz"
snp_tbi="${snp_vcf}.tbi"

download_if_missing "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz" "$snp_txt"
download_if_missing "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.sql" "$snp_sql"

if [[ ! -s "$snp_vcf" || ! -s "$snp_tbi" ]]; then
  echo "Converting UCSC snp151Common table to FACETS SNP VCF"
  {
    printf "##fileformat=VCFv4.2\n"
    printf "##source=UCSC_hg38_snp151Common_converted_for_FACETS\n"
    printf "##reference=Homo_sapiens_assembly38.fasta\n"
    printf "##INFO=<ID=UCSC_SNP151_COMMON,Number=0,Type=Flag,Description=\"Variant converted from UCSC hg38 snp151Common table for FACETS snp-pileup\">\n"
    awk 'BEGIN{FS="\t"} {print "##contig=<ID="$1",length="$2">"}' "$ref_fai"
    printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    gzip -dc "$snp_txt" | awk 'BEGIN{FS=OFS="\t"}
      function upper_base(x) { x=toupper(x); return x }
      {
        chrom=$2; start=$3; id=$5; ref=upper_base($9); observed=upper_base($10);
        class=$12; locType=$17; exceptions=$19;
        if (class != "single" || locType != "exact" || exceptions != "") next;
        if (ref !~ /^[ACGT]$/) next;
        n=split(observed, alleles, "/");
        if (n != 2) next;
        a1=upper_base(alleles[1]); a2=upper_base(alleles[2]);
        if (a1 !~ /^[ACGT]$/ || a2 !~ /^[ACGT]$/ || a1 == a2) next;
        if (a1 == ref) alt=a2; else if (a2 == ref) alt=a1; else next;
        print chrom, start + 1, id, ref, alt, ".", "PASS", "UCSC_SNP151_COMMON"
      }'
  } | bgzip -c > "${snp_vcf}.tmp"
  mv "${snp_vcf}.tmp" "$snp_vcf"
  tabix -f -p vcf "$snp_vcf"
fi

{
  printf "resource\tpath\n"
  printf "WES_SUBCLONE_REF_FASTA\t%s\n" "$ref_fasta"
  printf "WES_SUBCLONE_REF_CACHE\t%s\n" "${ref_cache_dir}/%2s/%2s/%s"
  printf "WES_SUBCLONE_REF_FASTA_FAI\t%s\n" "$ref_fai"
  printf "WES_SUBCLONE_REF_DICT\t%s\n" "$ref_dict"
  printf "WES_SUBCLONE_FACETS_SNP_VCF\t%s\n" "$snp_vcf"
  printf "WES_SUBCLONE_FACETS_SNP_VCF_TBI\t%s\n" "$snp_tbi"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${log_dir}/Auto_wes_subclone_resources.tsv"

echo "Prepared resources:"
cat "${log_dir}/Auto_wes_subclone_resources.tsv"
