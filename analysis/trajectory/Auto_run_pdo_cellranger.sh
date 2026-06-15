#!/bin/bash
#PBS -l select=1:ncpus=16:mem=512gb
#PBS -l walltime=48:00:00
#PBS -N Auto_PDO_VelCR
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<SUR_ID_PDO>"
  exit 1
fi

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
OUT="${WD}/PDOs_outs/Auto_velocity_PDO"
manifest="${OUT}/tables/Auto_pdo_velocity_sample_manifest.csv"
cellranger="/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-9.0.1/bin/cellranger"
transcriptome="/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A"

if [[ ! -f "$manifest" ]]; then
  echo "ERROR: missing manifest: $manifest"
  exit 1
fi

line=$(awk -F, -v s="$sample" 'NR > 1 && $1 == s {print; exit}' "$manifest")
if [[ -z "$line" ]]; then
  echo "ERROR: sample not found in manifest: $sample"
  exit 1
fi

batch_type=$(awk -F, '{print $2}' <<< "$line")
fastq_dir=$(awk -F, '{print $5}' <<< "$line")
run_root="${OUT}/cellranger"
run_dir="${run_root}/${sample}"

if [[ "$batch_type" != "Cynthia_batch" ]]; then
  echo "Sample $sample is not Cynthia_batch; CellRanger rerun is not required."
  exit 0
fi
if [[ ! -d "$fastq_dir" ]]; then
  echo "ERROR: missing FASTQ directory: $fastq_dir"
  exit 1
fi
if [[ ! -x "$cellranger" ]]; then
  chmod +x "$cellranger" || true
fi
if [[ ! -x "$cellranger" ]]; then
  echo "ERROR: CellRanger executable is not executable: $cellranger"
  exit 1
fi
if [[ -f "${run_dir}/outs/possorted_genome_bam.bam" ]]; then
  echo "Existing CellRanger BAM found for $sample; skipping CellRanger."
  echo $(date +%T)
  exit 0
fi
if [[ -e "$run_dir" ]]; then
  echo "ERROR: CellRanger directory exists but BAM is missing: $run_dir"
  echo "Refusing to overwrite or delete existing outputs."
  exit 1
fi

mkdir -p "$run_root" "${OUT}/logs"

sample_arg=$(
  find "$fastq_dir" -maxdepth 1 -type f -name "*.fastq.gz" -printf "%f\n" \
    | sed -E 's/_S[0-9]+_L[0-9]+_[IR][12]_001\.fastq\.gz$//' \
    | sort -u \
    | paste -sd, -
)

if [[ -z "$sample_arg" ]]; then
  echo "ERROR: no FASTQ sample prefix detected under $fastq_dir"
  exit 1
fi

echo "Sample: $sample"
echo "FASTQs: $fastq_dir"
echo "CellRanger sample argument: $sample_arg"
echo "Output root: $run_root"

cd "$run_root"
"$cellranger" count \
  --id="$sample" \
  --fastqs="$fastq_dir" \
  --transcriptome="$transcriptome" \
  --sample="$sample_arg" \
  --create-bam=true \
  --localcores=16 \
  --localmem=460

echo $(date +%T)
