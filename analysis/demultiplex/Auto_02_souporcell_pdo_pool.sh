#!/bin/bash
#PBS -l select=1:ncpus=18:mem=512gb
#PBS -l walltime=24:00:00
#PBS -N Auto_PDO_Souporcell
#PBS -koed

####################
# Analysis registry
# Status: active pool-level Souporcell PBS workflow.
# Script: analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh
# Methodology: analysis/methodology/demultiplex/demultiplex_methodology.md
# Map: analysis/ANALYSIS_MAP.md
# Inputs: Cell Ranger BAM, BAM index, filtered barcode TSV, matching genome
#         FASTA, and the live Souporcell/Demuxafy container. Historical BAM and
#         barcode paths can be supplied explicitly through PBS variables.
# Outputs: ephemeral demultiplex_intermediate/souporcell/<pool>/ for the full
#          heavy run; persistent live
#          souporcell_assignments/<pool>/Auto_<pool>_clusters.tsv.
# Downstream: cluster_genotypes.vcf feeds Auto_03_reference_and_assign.sh;
#             clusters.tsv feeds donor-specific count export.
####################

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

out_root="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/demultiplex_intermediate"
cellranger_out="${out_root}/cellranger/${pool}/outs"
souporcell_out="${out_root}/souporcell/${pool}"
genome_fasta="/rds/general/project/tumourheterogeneity1/live/demultiplex/genome.fa"
demuxafy_sif="/rds/general/project/spatialtranscriptomics/live/multiplexed/Demuxafy.sif"
souporcell_sif="/rds/general/project/tumourheterogeneity1/live/demultiplex/souporcell_latest.sif"
tmp_root="/rds/general/user/sg3723/home/tmpfiles"

####################
# Persist the small downstream-critical cell/cluster assignment table in live
# storage while retaining the full heavy Souporcell run under ephemeral.
live_demultiplex_root="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/demultiplex"
live_souporcell_dir="${live_demultiplex_root}/souporcell_assignments/${pool}"
live_clusters_tsv="${live_souporcell_dir}/Auto_${pool}_clusters.tsv"
####################

####################
# Use an explicitly supplied reference FASTA when a delivery requires the
# matching vendor reference; retain the historical reference by default.
input_genome_fasta="${input_genome_fasta:-}"
if [[ -n "$input_genome_fasta" ]]; then
  genome_fasta="$input_genome_fasta"
fi
genome_fasta_dir="$(dirname "$genome_fasta")"
####################

bam="${cellranger_out}/possorted_genome_bam.bam"
barcodes_gz="${cellranger_out}/filtered_feature_bc_matrix/barcodes.tsv.gz"
barcodes_tsv="${souporcell_out}/barcodes.tsv"

####################
# Permit reuse of retained historical Cell Ranger products without copying
# 20+ GB BAMs into the organized ephemeral rerun tree.
input_bam="${input_bam:-}"
input_barcodes_gz="${input_barcodes_gz:-}"
if [[ -n "$input_bam" ]]; then
  bam="$input_bam"
fi
if [[ -n "$input_barcodes_gz" ]]; then
  barcodes_gz="$input_barcodes_gz"
fi
input_bam_dir="$(dirname "$bam")"
input_barcodes_dir="$(dirname "$barcodes_gz")"
####################

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
  --bind "$input_bam_dir" \
  --bind "$input_barcodes_dir" \
  --bind "$genome_fasta_dir" \
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

####################
if [[ ! -s "${souporcell_out}/clusters.tsv" ]]; then
  echo "ERROR: Souporcell completed without clusters.tsv: ${souporcell_out}/clusters.tsv"
  exit 1
fi
mkdir -p "$live_souporcell_dir"
cp --preserve=timestamps "${souporcell_out}/clusters.tsv" "$live_clusters_tsv"
echo "Persistent cluster assignments: $live_clusters_tsv"
####################

echo $(date +%T)
