#!/bin/bash
#PBS -l select=1:ncpus=18:mem=512gb
#PBS -l walltime=24:00:00
#PBS -N Auto_PDO_Souporcell
#PBS -koed

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

pool="${pool:-}"
k="${k:-}"

if [[ -z "$pool" ]]; then
  echo "ERROR: submit with -v pool=PDOs_Untreated,k=6 or -v pool=PDOs_Treated,k=4"
  exit 1
fi
if [[ -z "$k" ]]; then
  case "$pool" in
    PDOs_Untreated) k=6 ;;
    PDOs_Treated) k=4 ;;
    *)
      echo "ERROR: k was not supplied and no default exists for pool=$pool"
      exit 1
      ;;
  esac
fi

out_root="/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex"
cellranger_out="${out_root}/cellranger/${pool}/outs"
souporcell_out="${out_root}/souporcell/${pool}"
genome_fasta="/rds/general/project/tumourheterogeneity1/live/demultiplex/genome.fa"
demuxafy_sif="/rds/general/project/spatialtranscriptomics/live/multiplexed/Demuxafy.sif"
souporcell_sif="/rds/general/project/tumourheterogeneity1/live/demultiplex/souporcell_latest.sif"
tmp_root="/rds/general/user/sg3723/home/tmpfiles"

bam="${cellranger_out}/possorted_genome_bam.bam"
barcodes_gz="${cellranger_out}/filtered_feature_bc_matrix/barcodes.tsv.gz"
barcodes_tsv="${souporcell_out}/barcodes.tsv"

if [[ ! -f "$bam" ]]; then
  echo "ERROR: missing CellRanger BAM: $bam"
  exit 1
fi
if [[ ! -f "$barcodes_gz" ]]; then
  echo "ERROR: missing filtered barcodes: $barcodes_gz"
  exit 1
fi
if [[ ! -f "$genome_fasta" ]]; then
  echo "ERROR: missing genome FASTA: $genome_fasta"
  exit 1
fi

container="$demuxafy_sif"
if [[ ! -f "$container" ]]; then
  container="$souporcell_sif"
fi
if [[ ! -f "$container" ]]; then
  echo "ERROR: no Souporcell/Demuxafy container found"
  exit 1
fi

mkdir -p "$souporcell_out" "$tmp_root" "${out_root}/logs"
zcat "$barcodes_gz" > "$barcodes_tsv"

export TMPDIR="$tmp_root"
echo "Pool: $pool"
echo "k: $k"
echo "BAM: $bam"
echo "Barcodes copy: $barcodes_tsv"
echo "Output: $souporcell_out"
echo "Container: $container"
command -v singularity

singularity exec \
  --bind "$out_root" \
  --bind "/rds/general/project/tumourheterogeneity1/live/demultiplex" \
  --bind "$tmp_root" \
  --bind /tmp \
  -B "$TMPDIR" \
  "$container" \
  souporcell_pipeline.py \
    -i "$bam" \
    -b "$barcodes_tsv" \
    -f "$genome_fasta" \
    -t 18 \
    -o "$souporcell_out" \
    -k "$k"

echo $(date +%T)
