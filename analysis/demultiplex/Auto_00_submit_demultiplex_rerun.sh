#!/bin/bash

set -euo pipefail

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
cd "$WD"

current_jobs=$(qstat | grep sg3723 | wc -l)
if [[ "$current_jobs" -gt 46 ]]; then
  echo "ERROR: current PBS job count is $current_jobs, above the configured limit."
  exit 1
fi

submit_log="analysis/demultiplex/Auto_demultiplex_rerun_jobs.tsv"
mkdir -p "$(dirname "$submit_log")"
printf "step\tpool\tjob_id\tdependency\n" > "$submit_log"

submit_job() {
  local step="$1"
  local pool="$2"
  local dependency="$3"
  local script="$4"
  local vars="$5"
  local job_id

  if [[ -n "$dependency" ]]; then
    job_id=$(qsub -W "depend=afterok:${dependency}" -v "$vars" "$script")
  else
    job_id=$(qsub -v "$vars" "$script")
  fi
  printf "%s\t%s\t%s\t%s\n" "$step" "$pool" "$job_id" "${dependency:-none}" | tee -a "$submit_log"
}

cell_untreated=$(qsub -v pool=PDOs_Untreated analysis/demultiplex/Auto_01_cellranger_pdo_pool.sh)
printf "cellranger\tPDOs_Untreated\t%s\tnone\n" "$cell_untreated" | tee -a "$submit_log"

cell_treated=$(qsub -v pool=PDOs_Treated analysis/demultiplex/Auto_01_cellranger_pdo_pool.sh)
printf "cellranger\tPDOs_Treated\t%s\tnone\n" "$cell_treated" | tee -a "$submit_log"

soup_untreated=$(qsub -W "depend=afterok:${cell_untreated}" -v pool=PDOs_Untreated,k=6 analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh)
printf "souporcell\tPDOs_Untreated\t%s\t%s\n" "$soup_untreated" "$cell_untreated" | tee -a "$submit_log"

soup_treated=$(qsub -W "depend=afterok:${cell_treated}" -v pool=PDOs_Treated,k=4 analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh)
printf "souporcell\tPDOs_Treated\t%s\t%s\n" "$soup_treated" "$cell_treated" | tee -a "$submit_log"

assign_untreated=$(qsub -W "depend=afterok:${soup_untreated}" -v pool=PDOs_Untreated analysis/demultiplex/Auto_03_reference_and_assign.sh)
printf "assign\tPDOs_Untreated\t%s\t%s\n" "$assign_untreated" "$soup_untreated" | tee -a "$submit_log"

assign_treated=$(qsub -W "depend=afterok:${soup_treated}" -v pool=PDOs_Treated analysis/demultiplex/Auto_03_reference_and_assign.sh)
printf "assign\tPDOs_Treated\t%s\t%s\n" "$assign_treated" "$soup_treated" | tee -a "$submit_log"

write_untreated=$(qsub -W "depend=afterok:${assign_untreated}" -v pool=PDOs_Untreated analysis/demultiplex/Auto_04_write_demultiplexed_counts.sh)
printf "write_counts\tPDOs_Untreated\t%s\t%s\n" "$write_untreated" "$assign_untreated" | tee -a "$submit_log"

write_treated=$(qsub -W "depend=afterok:${assign_treated}" -v pool=PDOs_Treated analysis/demultiplex/Auto_04_write_demultiplexed_counts.sh)
printf "write_counts\tPDOs_Treated\t%s\t%s\n" "$write_treated" "$assign_treated" | tee -a "$submit_log"

verify_job=$(qsub -W "depend=afterok:${write_untreated}:${write_treated}" analysis/demultiplex/Auto_05_verify_existing_assignments.sh)
printf "verify\tall\t%s\t%s:%s\n" "$verify_job" "$write_untreated" "$write_treated" | tee -a "$submit_log"

echo "Submitted demultiplex rerun. Job log: $submit_log"
