#!/bin/bash
#PBS -l select=1:ncpus=4:mem=24gb
#PBS -l walltime=36:00:00
#PBS -N Auto_PhyloWGS
#PBS -koed
####################
# Per-sample FACETS/PyClone -> PhyloWGS PBS wrapper.
#
# Builds auditable PhyloWGS SSM/CNV inputs, then runs PhyloWGS multi-chain
# MCMC if the ephemeral PhyloWGS environment has been prepared.
####################

set -euo pipefail

echo $(date +%T)

SCRIPT_DIR="/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/analysis/cnv/wes_subclone"
source "${SCRIPT_DIR}/Auto_wes_subclone_config.sh"

sample="${sample:-}"
if [[ -z "$sample" ]]; then
  echo "ERROR: submit with -v sample=<sample>" >&2
  exit 1
fi
if ! wes_subclone_is_high_confidence "$sample"; then
  echo "ERROR: ${sample} is not in the high-confidence WES subclone allow-list." >&2
  exit 1
fi

EPHEMERAL_ROOT="${EPHEMERAL_OUT_ROOT:-/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline/PDOs_outs/Auto_wes_subclone}"
PHYLOWGS_ENV="${PHYLOWGS_ENV:-${EPHEMERAL_ROOT}/phylowgs_env}"
PHYLOWGS_DIR="${PHYLOWGS_DIR:-${EPHEMERAL_ROOT}/tools/phylowgs}"
PHYLOWGS_NUM_CHAINS="${PHYLOWGS_NUM_CHAINS:-4}"
PHYLOWGS_BURNIN="${PHYLOWGS_BURNIN:-1000}"
PHYLOWGS_MCMC="${PHYLOWGS_MCMC:-2500}"
PHYLOWGS_RUN_SUFFIX="${PHYLOWGS_RUN_SUFFIX:-}"

module purge
module load tools/dev

eval "$(~/miniforge3/bin/conda shell.bash hook)"

conda activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
cd "$PROJECT_DIR"
Rscript "${SCRIPT_DIR}/Auto_make_phylowgs_inputs.R" "$sample" "$OUT_ROOT" "$EPHEMERAL_ROOT"

if [[ ! -x "${PHYLOWGS_ENV}/bin/python2" || ! -f "${PHYLOWGS_DIR}/multievolve.py" || ! -x "${PHYLOWGS_DIR}/mh.o" ]]; then
  echo "ERROR: PhyloWGS env/tool is not ready. Submit Auto_prepare_phylowgs_env.sh first." >&2
  exit 1
fi

conda activate "$PHYLOWGS_ENV"
export LD_LIBRARY_PATH="${PHYLOWGS_ENV}/lib:${LD_LIBRARY_PATH:-}"

INPUT_DIR="${OUT_ROOT}/tables/phylowgs/${sample}"
RUN_DIR="${EPHEMERAL_ROOT}/intermediate/phylowgs/${sample}${PHYLOWGS_RUN_SUFFIX}"
LIVE_REPORT_DIR="${OUT_ROOT}/reports/phylowgs/${sample}"
LIVE_LOG_DIR="${OUT_ROOT}/logs"
mkdir -p "$RUN_DIR" "$LIVE_REPORT_DIR" "$LIVE_LOG_DIR"

cp "${INPUT_DIR}/ssm_data.txt" "${RUN_DIR}/ssm_data.txt"
cp "${INPUT_DIR}/cnv_data.txt" "${RUN_DIR}/cnv_data.txt"

cd "$RUN_DIR"
python2 "${PHYLOWGS_DIR}/multievolve.py" \
  --num-chains "$PHYLOWGS_NUM_CHAINS" \
  --ssms ssm_data.txt \
  --cnvs cnv_data.txt \
  --output-dir chains \
  --burnin-samples "$PHYLOWGS_BURNIN" \
  --mcmc-samples "$PHYLOWGS_MCMC"

trees_zip="$(find "$RUN_DIR/chains" -name trees.zip -type f | sort | head -n 1 || true)"
if [[ -z "$trees_zip" ]]; then
  echo "ERROR: PhyloWGS finished without a trees.zip under ${RUN_DIR}/chains" >&2
  exit 1
fi

mkdir -p "${RUN_DIR}/json"
cd "${RUN_DIR}/json"
python2 "${PHYLOWGS_DIR}/write_results.py" \
  "$sample" \
  "$trees_zip" \
  "${sample}.summ.json.gz" \
  "${sample}.muts.json.gz" \
  "${sample}.mutass.zip" \
  --include-ssm-names

cp "$trees_zip" "${LIVE_REPORT_DIR}/${sample}.trees.zip"
cp "${RUN_DIR}/json/${sample}.summ.json.gz" "${LIVE_REPORT_DIR}/"
cp "${RUN_DIR}/json/${sample}.muts.json.gz" "${LIVE_REPORT_DIR}/"
cp "${RUN_DIR}/json/${sample}.mutass.zip" "${LIVE_REPORT_DIR}/"

{
  printf "sample\t%s\n" "$sample"
  printf "input_dir\t%s\n" "$INPUT_DIR"
  printf "run_dir\t%s\n" "$RUN_DIR"
  printf "run_suffix\t%s\n" "$PHYLOWGS_RUN_SUFFIX"
  printf "phylowgs_dir\t%s\n" "$PHYLOWGS_DIR"
  printf "num_chains\t%s\n" "$PHYLOWGS_NUM_CHAINS"
  printf "burnin_samples\t%s\n" "$PHYLOWGS_BURNIN"
  printf "mcmc_samples\t%s\n" "$PHYLOWGS_MCMC"
  printf "trees_zip\t%s\n" "${LIVE_REPORT_DIR}/${sample}.trees.zip"
  printf "summary_json\t%s\n" "${LIVE_REPORT_DIR}/${sample}.summ.json.gz"
  printf "mutations_json\t%s\n" "${LIVE_REPORT_DIR}/${sample}.muts.json.gz"
  printf "mutation_assignment_zip\t%s\n" "${LIVE_REPORT_DIR}/${sample}.mutass.zip"
  printf "finished\t%s\n" "$(date -Iseconds)"
} > "${LIVE_LOG_DIR}/Auto_${sample}_phylowgs_run.tsv"

echo $(date +%T)
