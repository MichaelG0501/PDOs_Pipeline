#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb
#PBS -l walltime=24:00:00
#PBS -N Auto_HATCHetCNA
#PBS -koed
####################
# Analysis registry
# Status: active
# Script: analysis/cnv/wes_subclone/Auto_run_hatchet_clone_cna_sample.sh
# Methodology: analysis/methodology/cnv/wes_subclone/Auto_wes_subclone_methodology.md
# Map: analysis/ANALYSIS_MAP.md CNV section.
# Inputs:
#   - CNVkit tumour .cnr target/antitarget bins from Sarek
#   - FACETS snp-pileup germline-heterozygous allele counts
#   - ephemeral HATCHet conda environment
# Outputs:
#   - live PDOs_outs/Auto_wes_clone_cna/tables/hatchet/<sample>/best.*.ucn and chosen/result UCN tables
#   - live PDOs_outs/Auto_wes_clone_cna/figures/hatchet/<sample>/
#   - live PDOs_outs/Auto_wes_clone_cna/logs/Auto_<sample>_hatchet_clone_cna_status.tsv
# Downstream use: model-based HATCHet/HATCHet2 WES clone-specific CNA inference.
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
  echo "ERROR: HATCHet clone-CNA run is restricted to high-confidence WES pairs." >&2
  exit 1
fi

LIVE_PROJECT="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline"
EPHEMERAL_PROJECT="/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline"
LIVE_OUT="${LIVE_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
EPHEMERAL_OUT="${EPHEMERAL_PROJECT}/PDOs_outs/Auto_wes_clone_cna"
HATCHET_ENV_DIR="${EPHEMERAL_OUT}/hatchet_env"
DMTCP_RSCRIPT="/rds/general/user/sg3723/home/anaconda3/envs/dmtcp/bin/Rscript"
BIN_SIZE_BP="${HATCHET_BIN_SIZE_BP:-250000}"
HATCHET_CLONE_RANGE="${HATCHET_CLONE_RANGE:-2,4}"
HATCHET_SEEDS="${HATCHET_SEEDS:-100}"
HATCHET_USE_FACETS_PURITY="${HATCHET_USE_FACETS_PURITY:-1}"
INTERMEDIATE_DIR="${EPHEMERAL_OUT}/intermediate/${sample}/hatchet_250kb"
LIVE_TABLE_DIR="${LIVE_OUT}/tables/hatchet/${sample}"
LIVE_FIG_DIR="${LIVE_OUT}/figures/hatchet/${sample}"
LOG_DIR="${LIVE_OUT}/logs"
INPUT_DIR="${LIVE_OUT}/tables/hatchet_inputs"
STATUS_FILE="${LOG_DIR}/Auto_${sample}_hatchet_clone_cna_status.tsv"
RUN_LOG="${LOG_DIR}/Auto_${sample}_hatchet_compute_cn.log"
CLUSTER_LOG="${LOG_DIR}/Auto_${sample}_hatchet_cluster_bins.log"
PLOT_LOG="${LOG_DIR}/Auto_${sample}_hatchet_plot_cn.log"

mkdir -p "$INTERMEDIATE_DIR" "$LIVE_TABLE_DIR" "$LIVE_FIG_DIR" "$LOG_DIR" "$INPUT_DIR"

write_status() {
  local status="$1"
  local message="$2"
  {
    printf "field\tvalue\n"
    printf "sample\t%s\n" "$sample"
    printf "status\t%s\n" "$status"
    printf "message\t%s\n" "$message"
    printf "hatchet_env_dir\t%s\n" "$HATCHET_ENV_DIR"
    printf "intermediate_dir\t%s\n" "$INTERMEDIATE_DIR"
    printf "live_table_dir\t%s\n" "$LIVE_TABLE_DIR"
    printf "live_fig_dir\t%s\n" "$LIVE_FIG_DIR"
    printf "bin_size_bp\t%s\n" "$BIN_SIZE_BP"
    printf "hatchet_clone_range\t%s\n" "$HATCHET_CLONE_RANGE"
    printf "hatchet_seeds\t%s\n" "$HATCHET_SEEDS"
    printf "hatchet_use_facets_purity\t%s\n" "$HATCHET_USE_FACETS_PURITY"
    printf "facets_purity\t%s\n" "${FACETS_PURITY:-NA}"
    printf "finished\t%s\n" "$(date -Iseconds)"
  } > "$STATUS_FILE"
}

if [[ ! -x "${HATCHET_ENV_DIR}/bin/hatchet" ]]; then
  write_status "blocked_missing_env" "Run Auto_prepare_hatchet_clone_cna_env.sh first; missing ${HATCHET_ENV_DIR}/bin/hatchet"
  cat "$STATUS_FILE"
  exit 1
fi
if [[ ! -x "$DMTCP_RSCRIPT" ]]; then
  write_status "blocked_missing_env" "Missing dmtcp Rscript: ${DMTCP_RSCRIPT}"
  cat "$STATUS_FILE"
  exit 1
fi

bb_file="${INPUT_DIR}/Auto_${sample}_hatchet_250kb.bb"
bb_summary="${LOG_DIR}/Auto_${sample}_hatchet_bb_summary.tsv"
"$DMTCP_RSCRIPT" "${SCRIPT_DIR}/Auto_make_hatchet_bb.R" "$sample" "$bb_file" "$bb_summary" "$BIN_SIZE_BP"

FACETS_PURITY="NA"
facets_purity_file="${LIVE_PROJECT}/PDOs_outs/Auto_wes_subclone/tables/facets/Auto_${sample}_facets_purity_ploidy.tsv"
if [[ -s "$facets_purity_file" ]]; then
  FACETS_PURITY="$(awk 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="purity") c=i; next} NR==2 && c{print $c}' "$facets_purity_file")"
fi
purity_args=()
if [[ "$HATCHET_USE_FACETS_PURITY" == "1" && "$FACETS_PURITY" != "NA" && -n "$FACETS_PURITY" ]]; then
  purity_args=(-P "$FACETS_PURITY")
fi

cp "$bb_file" "${INTERMEDIATE_DIR}/Auto_${sample}.bb"
cd "$INTERMEDIATE_DIR"

eval "$(~/miniforge3/bin/conda shell.bash hook)"
set +u
source activate "$HATCHET_ENV_DIR"
export LD_LIBRARY_PATH="${HATCHET_ENV_DIR}/lib:${LD_LIBRARY_PATH:-}"
####################
# HATCHet overlays a local hatchet.ini from the current working directory.
# The packaged C++ solver is linked to a missing Gurobi runtime on CX3, so use
# Pyomo/CBC from the same conda environment for compute-cn.
####################
cat > hatchet.ini <<'EOF'
[compute_cn]
solver = cbc
EOF
if [[ ! -x "${HATCHET_ENV_DIR}/bin/cbc" ]]; then
  write_status "blocked_missing_cbc" "Run Auto_prepare_hatchet_cbc_solver.sh first; missing ${HATCHET_ENV_DIR}/bin/cbc"
  cat "$STATUS_FILE"
  exit 1
fi
hatchet --help > "${LOG_DIR}/Auto_${sample}_hatchet_help.txt" 2>&1 || true

hatchet cluster-bins "Auto_${sample}.bb" \
  -o "Auto_${sample}.seg" \
  -O "Auto_${sample}.bbc" \
  -d 0.08 \
  --exactK 12 \
  -R 10 \
  --allow_gaps \
  > "$CLUSTER_LOG" 2>&1

hatchet compute-cn \
  -i "Auto_${sample}" \
  -n"$HATCHET_CLONE_RANGE" \
  "${purity_args[@]}" \
  -p "$HATCHET_SEEDS" \
  -j "${PBS_NP:-8}" \
  -v 2 \
  -u 0.06 \
  -r 12 \
  -eD 6 \
  -eT 12 \
  -tB 0.03 \
  -tR 0.5 \
  -l 0.5 \
  > "$RUN_LOG" 2>&1

if [[ -s best.bbc.ucn ]]; then
  hatchet plot-cn best.bbc.ucn -rC 10 -rG 1 > "$PLOT_LOG" 2>&1 || true
fi
conda deactivate
set -u

if [[ ! -s best.seg.ucn || ! -s best.bbc.ucn ]]; then
  write_status "failed_no_best_ucn" "HATCHet compute-cn completed without non-empty best.seg.ucn/best.bbc.ucn."
  cat "$STATUS_FILE"
  exit 1
fi

find "$INTERMEDIATE_DIR" -maxdepth 1 -type f \
  \( -name '*.ucn' -o -name 'results.*.tsv' -o -name 'Auto_*.bbc' -o -name 'Auto_*.seg' -o -name 'Auto_*.bb' \) \
  -exec cp {} "$LIVE_TABLE_DIR/" \;
find "$INTERMEDIATE_DIR" -maxdepth 1 -type f \
  \( -name '*.pdf' -o -name '*.png' \) \
  -exec cp {} "$LIVE_FIG_DIR/" \; || true

write_status "complete" "HATCHet completed and produced best.seg.ucn/best.bbc.ucn."
cat "$STATUS_FILE"
echo $(date +%T)
