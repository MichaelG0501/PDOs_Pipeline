#!/bin/bash
#PBS -l select=1:ncpus=16:mem=512gb
#PBS -l walltime=48:00:00
#PBS -N Auto_PDO_CellRanger
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

pool="${pool:-}"
if [[ -z "$pool" ]]; then
  echo "ERROR: submit with -v pool=PDOs_Untreated or -v pool=PDOs_Treated"
  exit 1
fi

raw_root="/rds/general/project/tumourheterogeneity1/live/ITH_sc/X204SC25083484-Z01-F001/X204SC25083484-Z01-F001/01.RawData"
out_root="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex"
cellranger="/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-9.0.1/bin/cellranger"
transcriptome="/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A"

fastq_dir="${raw_root}/${pool}"
run_root="${out_root}/cellranger"
run_dir="${run_root}/${pool}"

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
if [[ -e "$run_dir" ]]; then
  echo "ERROR: clean output directory already exists: $run_dir"
  echo "Refusing to overwrite or delete existing outputs."
  exit 1
fi

mkdir -p "$run_root" "${out_root}/logs"

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

echo "Pool: $pool"
echo "FASTQs: $fastq_dir"
echo "CellRanger sample argument: $sample_arg"
echo "Output root: $run_root"

cd "$run_root"
"$cellranger" count \
  --id="$pool" \
  --fastqs="$fastq_dir" \
  --transcriptome="$transcriptome" \
  --sample="$sample_arg" \
  --create-bam=true \
  --localcores=16 \
  --localmem=460

echo $(date +%T)
