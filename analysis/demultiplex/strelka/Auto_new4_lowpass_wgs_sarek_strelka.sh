#!/bin/bash
#PBS -l select=1:ncpus=32:mem=256gb
#PBS -l walltime=48:00:00
#PBS -N Auto_new4_Strelka
#PBS -koed

####################
# Analysis registry
# Status: active upstream donor-genotype generation workflow.
# Script: analysis/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_strelka.sh
# Methodology: analysis/methodology/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_strelka_methodology.md
# Map: analysis/ANALYSIS_MAP.md
# Inputs: analysis/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_samplesheet.csv
#         and its eight paired low-pass WGS FASTQ entries for SUR1346,
#         SUR1363, SUR1384, and SUR1391.
# Outputs: ephemeral Sarek/Nextflow results and work files under
#          PDOs_outs/demultiplex_intermediate/strelka_sarek/new4samples/;
#          persistent donor VCFs and indexes under the requested live
#          spatialtranscriptomics sarek_mutect/variant_calling/strelka/ tree;
#          a lightweight run manifest under live PDOs_outs/demultiplex/.
# Downstream: the four *.strelka.variants.vcf.gz files are reference donor
#             genotype inputs for Auto_03_reference_and_assign.sh.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/sarek

wd="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
samplesheet="${wd}/analysis/demultiplex/strelka/Auto_new4_lowpass_wgs_sarek_samplesheet.csv"
ephemeral_run_root="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/demultiplex_intermediate/strelka_sarek/new4samples"
sarek_out="${ephemeral_run_root}/results"
work_dir="${ephemeral_run_root}/work"
live_strelka_root="/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/strelka"
live_log_dir="${wd}/PDOs_outs/demultiplex/strelka/logs"
run_summary="${live_log_dir}/Auto_new4_lowpass_wgs_sarek_strelka_run_summary.tsv"
samples=(SUR1346 SUR1363 SUR1384 SUR1391)

if [[ ! -f "$samplesheet" ]]; then
  echo "ERROR: missing Sarek samplesheet: $samplesheet"
  exit 1
fi
if [[ ! -x /rds/general/user/sg3723/home/anaconda3/envs/sarek/bin/nextflow ]]; then
  echo "ERROR: nextflow is unavailable in the Sarek environment"
  exit 1
fi

while IFS=, read -r patient status sample lane fastq_1 fastq_2; do
  if [[ "$patient" == "patient" ]]; then
    continue
  fi
  if [[ "$status" != "0" ]]; then
    echo "ERROR: expected status=0 for germline-style Strelka calling: $sample $lane"
    exit 1
  fi
  if [[ ! -f "$fastq_1" || ! -f "$fastq_2" ]]; then
    echo "ERROR: missing FASTQ pair for $sample $lane"
    exit 1
  fi
done < "$samplesheet"

for sample in "${samples[@]}"; do
  destination_dir="${live_strelka_root}/${sample}"
  destination_vcf="${destination_dir}/${sample}.strelka.variants.vcf.gz"
  if [[ -e "$destination_vcf" || -e "${destination_vcf}.tbi" ]]; then
    echo "ERROR: refusing to overwrite existing donor VCF or index: $destination_vcf"
    exit 1
  fi
done

mkdir -p "$ephemeral_run_root" "$sarek_out" "$work_dir" "$live_log_dir"
cd "$ephemeral_run_root"

nextflow run nf-core/sarek \
  -r 3.4.4 \
  -profile singularity \
  -w "$work_dir" \
  -resume \
  --input "$samplesheet" \
  --outdir "$sarek_out" \
  --genome GATK.GRCh38 \
  --trim_fastq \
  --tools strelka \
  --max_cpus 32 \
  --max_memory 240.GB \
  --max_time 47.h

for sample in "${samples[@]}"; do
  source_dir="${sarek_out}/variant_calling/strelka/${sample}"
  source_vcf="${source_dir}/${sample}.strelka.variants.vcf.gz"
  source_index="${source_vcf}.tbi"
  destination_dir="${live_strelka_root}/${sample}"

  if [[ ! -s "$source_vcf" || ! -s "$source_index" ]]; then
    echo "ERROR: Sarek completed without the required Strelka VCF/index for $sample"
    exit 1
  fi

  mkdir -p "$destination_dir"
  cp --preserve=timestamps "$source_vcf" "$source_index" "$destination_dir/"
done

{
  printf "field\tvalue\n"
  printf "status\tcompleted\n"
  printf "completed_at\t%s\n" "$(date --iso-8601=seconds)"
  printf "sarek_version\t3.4.4\n"
  printf "tools\tstrelka\n"
  printf "genome\tGATK.GRCh38\n"
  printf "input_samplesheet\t%s\n" "$samplesheet"
  printf "ephemeral_results\t%s\n" "$sarek_out"
  printf "ephemeral_work\t%s\n" "$work_dir"
  printf "persistent_strelka_root\t%s\n" "$live_strelka_root"
  for sample in "${samples[@]}"; do
    printf "persistent_vcf_%s\t%s/%s/%s.strelka.variants.vcf.gz\n" \
      "$sample" "$live_strelka_root" "$sample" "$sample"
  done
} > "$run_summary"

echo $(date +%T)
