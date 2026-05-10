#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=6:00:00
#PBS -N Auto_PDO_GenotypeAssign
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/prod
module load BCFtools/1.22-GCC-14.2.0

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

wd="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
pool="${pool:-}"
donors="${donors:-}"

if [[ -z "$pool" ]]; then
  echo "ERROR: submit with -v pool=PDOs_Untreated,donors=1070,1090,1072,1121,1141,1181"
  exit 1
fi
if [[ -z "$donors" ]]; then
  case "$pool" in
    PDOs_Untreated) donors="1070,1090,1072,1121,1141,1181" ;;
    PDOs_Treated) donors="1070,1090,1072,1181" ;;
    *)
      echo "ERROR: donors was not supplied and no default exists for pool=$pool"
      exit 1
      ;;
  esac
fi

out_root="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex"
strelka_root="/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/strelka"
souporcell_out="${out_root}/souporcell/${pool}"
ref_out="${out_root}/reference_genotypes/${pool}"
assign_out="${out_root}/genotype_assignment/${pool}"
cluster_vcf="${souporcell_out}/cluster_genotypes.vcf"

if [[ ! -f "$cluster_vcf" ]]; then
  echo "ERROR: missing Souporcell cluster VCF: $cluster_vcf"
  exit 1
fi

mkdir -p "$ref_out" "$assign_out"

IFS=',' read -r -a donor_array <<< "$donors"
per_donor_vcfs=()

for donor in "${donor_array[@]}"; do
  donor="${donor//[[:space:]]/}"
  src="${strelka_root}/NT_${donor}/NT_${donor}.strelka.variants.vcf.gz"
  tmp="${ref_out}/Auto_SUR${donor}.snps.raw.vcf.gz"
  het_tmp="${ref_out}/Auto_SUR${donor}.het.raw.vcf.gz"
  out="${ref_out}/Auto_SUR${donor}.het.vcf.gz"
  sample_file="${ref_out}/Auto_SUR${donor}.sample_name.txt"

  if [[ ! -f "$src" ]]; then
    echo "ERROR: missing normal Strelka VCF for donor SUR${donor}: $src"
    exit 1
  fi

  printf "SUR%s\n" "$donor" > "$sample_file"
  bcftools view \
    -v snps \
    -m2 \
    -M2 \
    -f PASS \
    -i 'FORMAT/DP>=4 && FORMAT/GQ>=20' \
    "$src" \
    -Oz \
    -o "$tmp"
  tabix -f -p vcf "$tmp"
  bcftools view -g het "$tmp" -Oz -o "$het_tmp"
  tabix -f -p vcf "$het_tmp"
  bcftools reheader -s "$sample_file" -o "$out" "$het_tmp"
  tabix -f -p vcf "$out"
  per_donor_vcfs+=("$out")
done

reference_vcf="${ref_out}/Auto_${pool}_merged.het.vcf.gz"
bcftools merge \
  --missing-to-ref \
  -m none \
  -Oz \
  -o "$reference_vcf" \
  "${per_donor_vcfs[@]}"
tabix -f -p vcf "$reference_vcf"

ref_gt="${assign_out}/Auto_${pool}_reference_gt.tsv"
cluster_gt="${assign_out}/Auto_${pool}_cluster_gt.tsv"
ref_samples="${assign_out}/Auto_${pool}_reference_samples.txt"
cluster_samples="${assign_out}/Auto_${pool}_cluster_samples.txt"

bcftools query -l "$reference_vcf" > "$ref_samples"
bcftools query -l "$cluster_vcf" > "$cluster_samples"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$reference_vcf" > "$ref_gt"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' "$cluster_vcf" > "$cluster_gt"

cd "$wd"
Rscript analysis/demultiplex/Auto_03_genotyping_save_assign.R \
  --pool "$pool" \
  --ref_gt "$ref_gt" \
  --cluster_gt "$cluster_gt" \
  --ref_samples "$ref_samples" \
  --cluster_samples "$cluster_samples" \
  --outdir "$assign_out"

echo $(date +%T)
