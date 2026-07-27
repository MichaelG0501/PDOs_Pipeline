#!/bin/bash

####################
# Analysis registry
# Status: active generic multiplexed-PDO submission entry point.
# Script: analysis/demultiplex/Auto_00_submit_demultiplex_pool.sh
# Methodology: analysis/methodology/demultiplex/demultiplex_methodology.md
# Map: analysis/ANALYSIS_MAP.md
# Inputs: a paired-FASTQ directory, matching reference FASTA, Souporcell k,
#         optional WES/VCF donor references in reference mode.
# Outputs: ephemeral Cell Ranger/Souporcell intermediates; live submission
#          manifest, assignment audit, and donor or provisional-cluster CSVs.
# Downstream: exported count CSVs are inputs to the PDO count-matrix workflow;
#             temporary labels must be replaced after genotype assignment.
#
# Required: pool, fastq_dir, k.
# Optional genotype_mode=reference (default) runs WES/VCF-based assignment;
# genotype_mode=temporary exports singlet Souporcell clusters under an explicit
# temporary label prefix and deliberately does not assert donor identity.
#
# Example for the four newly sequenced PDOs:
# bash analysis/demultiplex/Auto_00_submit_demultiplex_pool.sh \
#   new4samples /rds/general/project/tumourheterogeneity1/live/ITH_sc/new4samples 4 temporary
####################

set -euo pipefail

pool="${1:-}"
fastq_dir="${2:-}"
k="${3:-}"
genotype_mode="${4:-reference}"
count_output_dir="${5:-}"

if [[ -z "$pool" || -z "$fastq_dir" || -z "$k" ]]; then
  echo "Usage: $0 <pool> <fastq_dir> <souporcell_k> [reference|temporary] [count_output_dir]"
  exit 1
fi
if [[ "$genotype_mode" != "reference" && "$genotype_mode" != "temporary" ]]; then
  echo "ERROR: genotype_mode must be reference or temporary"
  exit 1
fi
if [[ ! -d "$fastq_dir" ]]; then
  echo "ERROR: FASTQ directory does not exist: $fastq_dir"
  exit 1
fi

####################
# Preserve the historical pools' default output arrangement.  The new4samples
# run uses a separate, clearly named count-matrix directory beside its FASTQs.
if [[ -z "$count_output_dir" ]]; then
  count_output_dir="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/demultiplex/counts_csv/${pool}"
fi
temporary_cluster_prefix="TEMP_${pool}_SouporcellCluster"
live_root="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
submission_dir="${live_root}/PDOs_outs/demultiplex/submissions"
mkdir -p "$submission_dir"
submission_log="${submission_dir}/Auto_${pool}_submission.tsv"
printf "step\tpool\tjob_id\tdependency\tmode\n" > "$submission_log"
####################

wait_for_capacity() {
  while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do
    sleep 180
  done
}

wait_for_capacity
cell_job=$(qsub -v "pool=${pool},input_fastq_dir=${fastq_dir}" analysis/demultiplex/Auto_01_cellranger_pdo_pool.sh)
printf "cellranger\t%s\t%s\tnone\t%s\n" "$pool" "$cell_job" "$genotype_mode" | tee -a "$submission_log"

wait_for_capacity
soup_job=$(qsub -W "depend=afterok:${cell_job}" -v "pool=${pool},k=${k},input_genome_fasta=${fastq_dir}/genome.fa" analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh)
printf "souporcell\t%s\t%s\t%s\t%s\n" "$pool" "$soup_job" "$cell_job" "$genotype_mode" | tee -a "$submission_log"

if [[ "$genotype_mode" == "reference" ]]; then
  wait_for_capacity
  assign_job=$(qsub -W "depend=afterok:${soup_job}" -v "pool=${pool}" analysis/demultiplex/Auto_03_reference_and_assign.sh)
  printf "genotype_assignment\t%s\t%s\t%s\t%s\n" "$pool" "$assign_job" "$soup_job" "$genotype_mode" | tee -a "$submission_log"

  wait_for_capacity
  write_job=$(qsub -W "depend=afterok:${assign_job}" -v "pool=${pool},counts_out_dir=${count_output_dir}" analysis/demultiplex/Auto_04_write_demultiplexed_counts.sh)
  printf "write_counts\t%s\t%s\t%s\t%s\n" "$pool" "$write_job" "$assign_job" "$genotype_mode" | tee -a "$submission_log"
else
  wait_for_capacity
  write_job=$(qsub -W "depend=afterok:${soup_job}" -v "pool=${pool},temporary_cluster_prefix=${temporary_cluster_prefix},counts_out_dir=${count_output_dir}" analysis/demultiplex/Auto_04_write_demultiplexed_counts.sh)
  printf "write_temporary_cluster_counts\t%s\t%s\t%s\t%s\n" "$pool" "$write_job" "$soup_job" "$genotype_mode" | tee -a "$submission_log"
fi

echo "Submitted ${pool}. Durable submission log: ${submission_log}"
echo "Count matrix destination: ${count_output_dir}"
