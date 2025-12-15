#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=8:00:00
#PBS -N master

echo $(date +%T)

WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline

cd $WD

missing_samples=()
done_samples=()
no_cell_samples=()

for sample_folder in PDOs_outs/by_samples/*_PDO/; do
  #while [[ $(qstat | wc -l) -gt 46 ]]; do
    #sleep 180
  #done

  sample=$(basename "$sample_folder")

  rds_file="PDOs_outs/by_samples/$sample/${sample}_rank4_9_nrun10.RDS"
  no_cell="PDOs_outs/by_samples/$sample/no_cell"

  if [[ ! -f "$rds_file" && ! -f "$no_cell" ]]; then
    qsub -v sample="$sample" -N "$sample" 2_NMF.sh
    missing_samples+=("$sample")
  else
    [[ -f "$rds_file" ]] && done_samples+=("$sample")
    [[ -f "$no_cell" ]] && no_cell_samples+=("$sample")
  fi
done

echo
echo "Jobs submitted (with epithelial cells): ${#missing_samples[@]}"
((${#missing_samples[@]})) && printf '  %s\n' "${missing_samples[@]}"

echo
echo "Completed (has RDS): ${#done_samples[@]}"
((${#done_samples[@]})) && printf '  %s\n' "${done_samples[@]}"

echo
echo "Skip marker no_cell: ${#no_cell_samples[@]}"
((${#no_cell_samples[@]})) && printf '  %s\n' "${no_cell_samples[@]}"

echo
echo "$(date +%T)"
