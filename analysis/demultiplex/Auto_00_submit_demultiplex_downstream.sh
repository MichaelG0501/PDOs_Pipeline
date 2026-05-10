#!/bin/bash

set -euo pipefail

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
cd "$WD"

submit_log="analysis/demultiplex/Auto_demultiplex_downstream_jobs.tsv"
printf "step\tpool\tjob_id\tdependency\n" > "$submit_log"

soup_untreated=$(qsub -v pool=PDOs_Untreated,k=6 analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh)
printf "souporcell\tPDOs_Untreated\t%s\tnone\n" "$soup_untreated" | tee -a "$submit_log"

soup_treated=$(qsub -v pool=PDOs_Treated,k=4 analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh)
printf "souporcell\tPDOs_Treated\t%s\tnone\n" "$soup_treated" | tee -a "$submit_log"

assign_untreated=$(qsub -W "depend=afterok:${soup_untreated}" -v pool=PDOs_Untreated analysis/demultiplex/Auto_03_reference_and_assign.sh)
printf "assign\tPDOs_Untreated\t%s\t%s\n" "$assign_untreated" "$soup_untreated" | tee -a "$submit_log"

assign_treated=$(qsub -W "depend=afterok:${soup_treated}" -v pool=PDOs_Treated analysis/demultiplex/Auto_03_reference_and_assign.sh)
printf "assign\tPDOs_Treated\t%s\t%s\n" "$assign_treated" "$soup_treated" | tee -a "$submit_log"

write_untreated=$(qsub -W "depend=afterok:${assign_untreated}" -v pool=PDOs_Untreated analysis/demultiplex/Auto_04_write_demultiplexed_counts.sh)
printf "write_counts\tPDOs_Untreated\t%s\t%s\n" "$write_untreated" "$assign_untreated" | tee -a "$submit_log"

write_treated=$(qsub -W "depend=afterok:${assign_treated}" -v pool=PDOs_Treated analysis/demultiplex/Auto_04_write_demultiplexed_counts.sh)
printf "write_counts\tPDOs_Treated\t%s\t%s\n" "$write_treated" "$assign_treated" | tee -a "$submit_log"

barcode_check=$(qsub -W "depend=afterok:${assign_untreated}" analysis/demultiplex/Auto_06_current_barcode_genotype_check.sh)
printf "current_barcode_genotype_check\tPDOs_Untreated\t%s\t%s\n" "$barcode_check" "$assign_untreated" | tee -a "$submit_log"

verify_job=$(qsub -W "depend=afterok:${write_untreated}:${write_treated}" analysis/demultiplex/Auto_05_verify_existing_assignments.sh)
printf "verify\tall\t%s\t%s:%s\n" "$verify_job" "$write_untreated" "$write_treated" | tee -a "$submit_log"

echo "Submitted downstream demultiplex rerun. Job log: $submit_log"
