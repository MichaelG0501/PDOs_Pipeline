#!/bin/bash
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -l walltime=02:00:00
#PBS -N Auto_ThetaN3
#PBS -koed
####################
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_import_theta2_n3_clone_cna.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - Sarek CNVkit tumour .cns for a high-confidence WES pair
#   - PDOs_outs/Auto_wes_clone_cna/tables/theta2/Auto_<sample>.n3.results
#   - ephemeral CNVkit Python 3.10 conda environment
# Outputs:
#   - live PDOs_outs/Auto_wes_clone_cna/tables/cns_theta2_n3/Auto_<sample>_theta2_n3_clone*.cns
#   - live PDOs_outs/Auto_wes_clone_cna/tables/Auto_theta2_n3_clone_cna_manifest.csv
#   - live PDOs_outs/Auto_wes_clone_cna/logs/Auto_<sample>_theta2_n3_import_status.tsv
# Downstream use: THetA2 forced n=3 clone-specific WES CNA comparison.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

####################
# Avoid inherited interactive conda state from qsub -V.
####################
unset CONDA_DEFAULT_ENV CONDA_PREFIX CONDA_PREFIX_1 CONDA_PROMPT_MODIFIER || true
unset CONDA_SHLVL CONDA_EXE CONDA_PYTHON_EXE _CE_CONDA _CE_M || true

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

sample="${sample:-PDO_1090_vs_NT_1090}"
if ! wes_subclone_is_high_confidence "$sample"; then
  echo "ERROR: THetA2 n3 import is restricted to high-confidence WES pairs." >&2
  exit 1
fi

LIVE_PROJECT="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
EPHEMERAL_PROJECT="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
LIVE_OUT="${LIVE_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
EPHEMERAL_OUT="${EPHEMERAL_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
CNVKIT_ENV_DIR="${EPHEMERAL_OUT}/cnvkit_py310_env"
CNS_DIR="${LIVE_OUT}/tables/cns_theta2_n3"
TABLE_DIR="${LIVE_OUT}/tables"
LOG_DIR="${LIVE_OUT}/logs"
IMPORT_DIR="${EPHEMERAL_OUT}/intermediate/${sample}/theta2_n3_import_${PBS_JOBID:-manual_$(date +%Y%m%d%H%M%S)}"
STATUS_FILE="${LOG_DIR}/Auto_${sample}_theta2_n3_import_status.tsv"
SAMPLE_MANIFEST="${LOG_DIR}/Auto_${sample}_theta2_n3_import_manifest.tsv"
GLOBAL_MANIFEST="${TABLE_DIR}/Auto_theta2_n3_clone_cna_manifest.csv"

mkdir -p "$CNS_DIR" "$TABLE_DIR" "$LOG_DIR" "$IMPORT_DIR"

write_status() {
  local status="$1"
  local message="$2"
  {
    printf "field\tvalue\n"
    printf "sample\t%s\n" "$sample"
    printf "status\t%s\n" "$status"
    printf "message\t%s\n" "$message"
    printf "cnvkit_env_dir\t%s\n" "$CNVKIT_ENV_DIR"
    printf "import_dir\t%s\n" "$IMPORT_DIR"
    printf "output_cns_dir\t%s\n" "$CNS_DIR"
    printf "finished\t%s\n" "$(date -Iseconds)"
  } > "$STATUS_FILE"
}

if [[ ! -x "${CNVKIT_ENV_DIR}/bin/cnvkit.py" ]]; then
  write_status "blocked_missing_env" "Missing CNVkit env: ${CNVKIT_ENV_DIR}/bin/cnvkit.py"
  cat "$STATUS_FILE"
  exit 1
fi

read -r tumor normal tumor_vcf_sample normal_vcf_sample tumor_cram normal_cram mutect_vcf cnvkit_cns < <(wes_subclone_paths_for_sample "$sample")
theta_results="${LIVE_OUT}/tables/theta2/Auto_${sample}.n3.results"
conditional_cns="${LIVE_PROJECT}/PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_${sample}_bulk_highres_conditional_shift_absolute.cns"
theta_cns="$cnvkit_cns"
if [[ -s "$conditional_cns" ]]; then
  theta_cns="$conditional_cns"
fi

for path in "$theta_cns" "$theta_results"; do
  if [[ ! -s "$path" ]]; then
    write_status "blocked_missing_input" "Missing required input: ${path}"
    cat "$STATUS_FILE"
    exit 1
  fi
done

eval "$(~/miniforge3/bin/conda shell.bash hook)"
set +u
source activate "$CNVKIT_ENV_DIR"
cnvkit.py import-theta "$theta_cns" "$theta_results" -d "$IMPORT_DIR"
conda deactivate
set -u

imported_files=( "$IMPORT_DIR"/*.cns )
if [[ ! -e "${imported_files[0]}" ]]; then
  write_status "failed_no_imported_cns" "CNVkit import-theta did not produce *.cns in ${IMPORT_DIR}"
  cat "$STATUS_FILE"
  exit 1
fi

{
  printf "sample\tclone_index\tsource_theta_results\timported_cns\toutput_cns\n"
  clone_index=0
  for imported in "${imported_files[@]}"; do
    clone_index=$((clone_index + 1))
    out_cns="${CNS_DIR}/Auto_${sample}_theta2_n3_clone${clone_index}_$(basename "$imported")"
    cp "$imported" "$out_cns"
    printf "%s\t%s\t%s\t%s\t%s\n" "$sample" "$clone_index" "$theta_results" "$imported" "$out_cns"
  done
} > "$SAMPLE_MANIFEST"

python - "$LOG_DIR" "$TABLE_DIR" "$GLOBAL_MANIFEST" <<'PY'
import csv
import glob
import os
import sys

log_dir, table_dir, global_manifest = sys.argv[1:4]
rows = []
for path in sorted(glob.glob(os.path.join(log_dir, "Auto_*_theta2_n3_import_manifest.tsv"))):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows.extend(reader)
os.makedirs(table_dir, exist_ok=True)
with open(global_manifest, "w", newline="") as handle:
    fieldnames = ["sample", "clone_index", "source_theta_results", "imported_cns", "output_cns"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
PY

write_status "complete" "Imported ${clone_index} THetA2 n3 clone CNS file(s)."
printf "theta_cns_input\t%s\n" "$theta_cns" >> "$STATUS_FILE"
cat "$STATUS_FILE"
echo $(date +%T)
