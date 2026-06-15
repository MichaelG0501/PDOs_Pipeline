#!/bin/bash
set -euo pipefail

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
OUT="${WD}/PDOs_outs/Auto_velocity_PDO"
cd "$WD"

mkdir -p "${OUT}/logs"

echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/trajectory/Auto_export_pdo_velocity_metadata.R

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate velocity
export PYTHONPYCACHEPREFIX="${OUT}/tmp_pycache"
python analysis/trajectory/Auto_prepare_pdo_velocity_refs.py

manifest="${OUT}/tables/Auto_pdo_velocity_sample_manifest.csv"
if [[ ! -f "$manifest" ]]; then
  echo "ERROR: missing manifest after metadata export: $manifest"
  exit 1
fi

throttle() {
  while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do
    sleep 180
  done
}

sanitize_job_name() {
  local x="$1"
  x="${x//[^A-Za-z0-9_]/_}"
  echo "${x:0:14}"
}

velocyto_jobs=()

while IFS=, read -r sample batch_type treatment pool fastq_dir cellranger_out bam barcodes_file n_cells has_bam; do
  [[ "$sample" == "sample" ]] && continue
  short_name=$(sanitize_job_name "$sample")

  if [[ "$batch_type" == "Cynthia_batch" ]]; then
    throttle
    jid_cr=$(qsub \
      -v sample="$sample" \
      -N "CR_${short_name}" \
      -o "${OUT}/logs/Auto_cellranger_${sample}.log" \
      -e "${OUT}/logs/Auto_cellranger_${sample}.err" \
      analysis/trajectory/Auto_run_pdo_cellranger.sh)
    filter_dep="-W depend=afterok:${jid_cr}"
    echo "Submitted CellRanger ${sample}: ${jid_cr}"
  else
    filter_dep=""
  fi

  throttle
  if [[ -n "$filter_dep" ]]; then
    jid_sort=$(qsub \
      $filter_dep \
      -v sample="$sample" \
      -N "Sort_${short_name}" \
      -o "${OUT}/logs/Auto_filter_sort_${sample}.log" \
      -e "${OUT}/logs/Auto_filter_sort_${sample}.err" \
      analysis/trajectory/Auto_filter_sort_pdo_velocity.sh)
  else
    jid_sort=$(qsub \
      -v sample="$sample" \
      -N "Sort_${short_name}" \
      -o "${OUT}/logs/Auto_filter_sort_${sample}.log" \
      -e "${OUT}/logs/Auto_filter_sort_${sample}.err" \
      analysis/trajectory/Auto_filter_sort_pdo_velocity.sh)
  fi
  echo "Submitted filter/sort ${sample}: ${jid_sort}"

  throttle
  jid_vel=$(qsub \
    -W depend=afterok:${jid_sort} \
    -v sample="$sample" \
    -N "Vel_${short_name}" \
    -o "${OUT}/logs/Auto_velocyto_${sample}.log" \
    -e "${OUT}/logs/Auto_velocyto_${sample}.err" \
    analysis/trajectory/Auto_run_pdo_velocyto.sh)
  echo "Submitted velocyto ${sample}: ${jid_vel}"
  velocyto_jobs+=("$jid_vel")
done < "$manifest"

dep=$(IFS=:; echo "${velocyto_jobs[*]}")
throttle
jid_vis=$(qsub \
  -W depend=afterok:${dep} \
  -o "${OUT}/logs/Auto_scvelo_visualisation.log" \
  -e "${OUT}/logs/Auto_scvelo_visualisation.err" \
  analysis/trajectory/Auto_run_pdo_scvelo_visualisation.sh)

echo "Submitted dependent scVelo visualisation: ${jid_vis}"
echo $(date +%T)
