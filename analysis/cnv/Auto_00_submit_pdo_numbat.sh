#!/bin/bashv/Auto_00_submit_pdo_numbat.sh
set -euo pipefail

WD="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
OUT="${WD}/PDOs_outs/Auto_PDO_numbat"
cd "$WD"

mkdir -p "${OUT}/logs"
run_tag=$(date +%Y%m%d_%H%M%S)
LOG_DIR="${OUT}/logs/Auto_numbat_run_${run_tag}"
mkdir -p "$LOG_DIR"

echo $(date +%T)
module purge
module load tools/dev

manifest="${OUT}/Auto_PDO_numbat_manifest.csv"
if [[ ! -f "$manifest" ]]; then
  eval "$(~/miniforge3/bin/conda shell.bash hook)"
  source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
  Rscript analysis/cnv/Auto_PDO_numbat_export_inputs.R
fi

if [[ ! -f "$manifest" ]]; then
  echo "ERROR: missing manifest after export: $manifest"
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

throttle
jid_img=$(qsub \
  -o "${LOG_DIR}/Auto_prepare_numbat_container.log" \
  -e "${LOG_DIR}/Auto_prepare_numbat_container.err" \
  analysis/cnv/Auto_prepare_pdo_numbat_container.sh)
echo "Submitted Numbat container preparation: ${jid_img}"
:
run_jobs=()
while IFS=, read -r sample batch_type treatment pool bam barcodes_file count_rds metadata_csv sample_out_dir numbat_dir allele_file clone_post_file joint_post_file n_cells has_
bam has_barcode_file; do
  [[ "$sample" == "sample" ]] && continue
  short_name=$(sanitize_job_name "$sample")
:
  throttle
  jid_pile=$(qsub \
    -W depend=afterok:${jid_img} \
    -v sample="$sample" \
    -N "NBp_${short_name}" \
    -o "${LOG_DIR}/Auto_numbat_pileup_${sample}.log" \
    -e "${LOG_DIR}/Auto_numbat_pileup_${sample}.err" \
    analysis/cnv/Auto_run_pdo_numbat_pileup.sh)
  echo "Submitted Numbat pileup ${sample}: ${jid_pile}"
:
  throttle
  jid_run=$(qsub \
    -W depend=afterok:${jid_pile} \
    -v sample="$sample" \
    -N "NBr_${short_name}" \
    -o "${LOG_DIR}/Auto_numbat_run_${sample}.log" \
    -e "${LOG_DIR}/Auto_numbat_run_${sample}.err" \
    analysis/cnv/Auto_run_pdo_numbat_sample.sh)
  echo "Submitted Numbat run ${sample}: ${jid_run}"
  run_jobs+=("$jid_run")
done < "$manifest"
:
dep=$(IFS=:; echo "${run_jobs[*]}")
throttle
jid_plot=$(qsub \
  -W depend=afterok:${dep} \
  -o "${LOG_DIR}/Auto_numbat_concordance_heatmaps.log" \
  -e "${LOG_DIR}/Auto_numbat_concordance_heatmaps.err" \
  analysis/cnv/Auto_run_pdo_numbat_concordance.sh)
echo "Submitted dependent concordance heatmaps: ${jid_plot}"
echo $(date +%T)
