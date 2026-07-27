#!/bin/bash
#PBS -l select=1:ncpus=8:mem=48gb
#PBS -l walltime=24:00:00
#PBS -N Auto_ThetaCNA
#PBS -koed
####################
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_run_theta2_clone_cna_sample.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - Sarek CNVkit tumour .cns and reference.cnn
#   - FACETS snp-pileup from the WES subclone workflow
#   - ephemeral THetA2 and CNVkit conda environments
# Outputs:
#   - live PDOs_outs/Auto_wes_clone_cna/tables/cns/Auto_<sample>_theta2*.cns
#   - live PDOs_outs/Auto_wes_clone_cna/tables/theta2/Auto_<sample>_theta2*.results
#   - live PDOs_outs/Auto_wes_clone_cna/logs/Auto_<sample>_theta2_clone_cna_status.tsv
# Downstream use: true model-based WES clone-specific CNA profiles when THetA2 succeeds.
####################

set -euo pipefail

echo $(date +%T)
module purge
module load tools/dev

####################
# Avoid inherited interactive conda state from qsub -V; those deactivate hooks
# can fail under set -u before the intended workflow environments activate.
####################
unset CONDA_DEFAULT_ENV CONDA_PREFIX CONDA_PREFIX_1 CONDA_PROMPT_MODIFIER || true
unset CONDA_SHLVL CONDA_EXE CONDA_PYTHON_EXE _CE_CONDA _CE_M || true

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

sample="${sample:-PDO_1181_vs_NT_1181}"
if ! wes_subclone_is_high_confidence "$sample"; then
  echo "ERROR: THetA2 clone-CNA run is currently restricted to high-confidence WES pairs." >&2
  exit 1
fi

LIVE_PROJECT="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
EPHEMERAL_PROJECT="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
LIVE_OUT="${LIVE_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
EPHEMERAL_OUT="${EPHEMERAL_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
THETA_ENV_DIR="${EPHEMERAL_OUT}/theta2_env"
CNVKIT_ENV_DIR="${EPHEMERAL_OUT}/cnvkit_py310_env"
DMTCP_RSCRIPT="/rds/general/user/sg3723/home/anaconda3/envs/dmtcp/bin/Rscript"
INTERMEDIATE_DIR="${EPHEMERAL_OUT}/intermediate/${sample}"
THETA_TABLE_DIR="${LIVE_OUT}/tables/theta2"
CNS_DIR="${LIVE_OUT}/tables/cns"
LOG_DIR="${LIVE_OUT}/logs"
STATUS_FILE="${LOG_DIR}/Auto_${sample}_theta2_clone_cna_status.tsv"
THETA_IMPORT_DIR="${INTERMEDIATE_DIR}/theta2_best_import_${PBS_JOBID:-manual_$(date +%Y%m%d%H%M%S)}"

mkdir -p "$INTERMEDIATE_DIR" "$THETA_TABLE_DIR" "$CNS_DIR" "$LOG_DIR" "$THETA_IMPORT_DIR"

write_status() {
  local status="$1"
  local message="$2"
  {
    printf "field\tvalue\n"
    printf "sample\t%s\n" "$sample"
    printf "status\t%s\n" "$status"
    printf "message\t%s\n" "$message"
    printf "theta_env_dir\t%s\n" "$THETA_ENV_DIR"
    printf "cnvkit_env_dir\t%s\n" "$CNVKIT_ENV_DIR"
    printf "intermediate_dir\t%s\n" "$INTERMEDIATE_DIR"
    printf "theta_import_dir\t%s\n" "$THETA_IMPORT_DIR"
    printf "finished\t%s\n" "$(date -Iseconds)"
  } > "$STATUS_FILE"
}

if [[ ! -x "${CNVKIT_ENV_DIR}/bin/cnvkit.py" || ! -x "${THETA_ENV_DIR}/bin/RunTHetA.py" ]]; then
  write_status "blocked_missing_env" "Run Auto_prepare_theta2_clone_cna_env.sh first; cnvkit.py or RunTHetA.py is unavailable."
  cat "$STATUS_FILE"
  exit 1
fi
if [[ ! -x "$DMTCP_RSCRIPT" ]]; then
  write_status "blocked_missing_env" "dmtcp Rscript is unavailable: ${DMTCP_RSCRIPT}"
  cat "$STATUS_FILE"
  exit 1
fi

eval "$(~/miniforge3/bin/conda shell.bash hook)"

read -r tumor normal tumor_vcf_sample normal_vcf_sample tumor_cram normal_cram mutect_vcf cnvkit_cns < <(wes_subclone_paths_for_sample "$sample")
cnvkit_dir="${SAREK_ROOT}/variant_calling/cnvkit/${sample}"
reference_cnn="${cnvkit_dir}/reference.cnn"
facets_pileup="${EPHEMERAL_PROJECT}/PDOs_outs/Auto_wes_subclone/intermediate/facets/${sample}/${sample}.snp_pileup.gz"
conditional_cns="${LIVE_PROJECT}/PDOs_outs/Auto_wes_absolute_cna/tables/cns/Auto_${sample}_bulk_highres_conditional_shift_absolute.cns"
theta_cns="$cnvkit_cns"
if [[ -s "$conditional_cns" ]]; then
  theta_cns="$conditional_cns"
fi

for path in "$theta_cns" "$reference_cnn" "$facets_pileup"; do
  if [[ ! -s "$path" ]]; then
    write_status "blocked_missing_input" "Missing required input: ${path}"
    cat "$STATUS_FILE"
    exit 1
  fi
done

interval_count="${INTERMEDIATE_DIR}/Auto_${sample}.interval_count"
normal_snp="${INTERMEDIATE_DIR}/Auto_${sample}.normal.snp_formatted.txt"
tumor_snp="${INTERMEDIATE_DIR}/Auto_${sample}.tumor.snp_formatted.txt"
snp_summary="${LOG_DIR}/Auto_${sample}_theta2_snp_counts.tsv"

if [[ "${PDO_FORCE_REBUILD:-0}" == "1" || ! -s "$interval_count" ]]; then
  set +u
  source activate "$CNVKIT_ENV_DIR"
  cnvkit.py export theta "$theta_cns" -r "$reference_cnn" -o "$interval_count"
  conda deactivate
  set -u
fi

if [[ "${PDO_FORCE_REBUILD:-0}" == "1" || ! -s "$normal_snp" || ! -s "$tumor_snp" ]]; then
  "$DMTCP_RSCRIPT" "${SCRIPT_DIR}/Auto_make_theta2_snp_counts.R" "$sample" "$facets_pileup" "$INTERMEDIATE_DIR" "$snp_summary"
fi

cd "$INTERMEDIATE_DIR"
set +u
source activate "$THETA_ENV_DIR"
RunTHetA.py "$interval_count" \
  --TUMOR_FILE "$tumor_snp" \
  --NORMAL_FILE "$normal_snp" \
  --BAF \
  --NUM_PROCESSES "${PBS_NP:-8}" \
  --FORCE
conda deactivate
set -u

best_results="$(find "$INTERMEDIATE_DIR" -maxdepth 1 -type f -name '*.BEST.results' -print | sort | tail -n 1)"
if [[ -z "$best_results" || ! -s "$best_results" ]]; then
  write_status "failed_no_best_results" "THetA2 completed without a non-empty BEST.results file."
  cat "$STATUS_FILE"
  exit 1
fi

cp "$best_results" "${THETA_TABLE_DIR}/Auto_${sample}_theta2.BEST.results"
find "$INTERMEDIATE_DIR" -maxdepth 1 -type f -name '*.results' -exec cp {} "$THETA_TABLE_DIR/" \;

set +u
source activate "$CNVKIT_ENV_DIR"
cnvkit.py import-theta "$theta_cns" "${THETA_TABLE_DIR}/Auto_${sample}_theta2.BEST.results" -d "$THETA_IMPORT_DIR"
conda deactivate
set -u

imported_files=( "$THETA_IMPORT_DIR"/*.cns )
if [[ ! -e "${imported_files[0]}" ]]; then
  write_status "failed_no_imported_cns" "CNVkit import-theta did not produce *.cns in ${THETA_IMPORT_DIR}"
  cat "$STATUS_FILE"
  exit 1
fi
imported_count=0
for cns in "${imported_files[@]}"; do
  imported_count=$((imported_count + 1))
  cp "$cns" "${CNS_DIR}/Auto_${sample}_theta2_best_clone${imported_count}_$(basename "$cns")"
done

write_status "complete" "THetA2 completed; imported ${imported_count} CNVkit .cns file(s)."
printf "theta_cns_input\t%s\n" "$theta_cns" >> "$STATUS_FILE"
cat "$STATUS_FILE"
echo $(date +%T)
