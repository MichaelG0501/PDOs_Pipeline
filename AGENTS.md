# AGENTS.md — PDOs_Pipeline

Single-cell RNA-seq QC and analysis of malignant OAC patient-derived organoids
on Imperial College HPC (PBS Pro). Core computation is R with bash PBS wrappers.
No cell-type annotation or CNA inference is required for the canonical PDO
state analysis.

## Authoritative documentation

- `analysis/ANALYSIS_MAP.md`: current run order, dependencies, status, and
  active/legacy classification.
- Script registry headers: exact file-level inputs, outputs, cache behavior,
  environment, and run command.
- `analysis/methodology/`: scientific rationale and thresholds for complex
  workflows.
- This file: durable repository-wide rules and canonical objects only.

When documents conflict, the current centred route in `ANALYSIS_MAP.md` and
the top authoritative script registry take precedence over historical comments.

## Core pipeline

| Step | Wrapper | R script | Scope | Environment |
| :--- | :--- | :--- | :--- | :--- |
| QC | `1_QC_Pipeline.sh` | `QC_Pipeline.R` | all samples | dmtcp |
| classical NMF (optional) | `2_master.sh` → `2_NMF.sh` | `NMF.R` | per sample | dmtcp |
| original GeneNMF | `3_geneNMF.sh` | `geneNMF.R` | all samples | gnmf |

The current downstream MP/state route is the centred workflow under
`analysis/metaprograms/centred/`, not the original uncentred nMP=13 route.
Always exclude `SUR843T3_PDO` from NMF and downstream analysis.

## Mandatory HPC and file-safety rules

1. Never run analytical, memory-heavy, or IO-heavy work on a login node. Use PBS.
   Login-node work is limited to small text/file inspections and syntax checks.
2. Every PBS job must include resource headers and `#PBS -koed` for live logs.
3. Initialize conda with
   `eval "$(~/miniforge3/bin/conda shell.bash hook)"` before activation.
4. Use `dmtcp` for general Seurat/analysis and `gnmf` for GeneNMF/UCell.
5. Keep at most 46 concurrent jobs; submitters throttle with
   `while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do sleep 180; done`.
6. All new persistent files start with `Auto_`. Superseded scripts use
   `legacy_`; manual-removal candidates use `delete_`. Do not delete files
   automatically.
7. Preserve user changes in dirty worktrees. Do not rewrite unrelated code.
8. New code added to an existing script is enclosed by 20-hash blocks:
   `####################`. Do not alter existing code outside such blocks
   without explicit permission.
9. Temporary R checks are named `delete_<description>.R` and removed
   immediately after use.
10. Analytical reruns requested by the user must be submitted and monitored
    until success and required outputs are verified.
11. Strictly forbid fallbacks: must use only the first and best option for any analysis. If not available or if it fails, the script should stop immediately.

PBS template:

```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
#PBS -koed
set -euo pipefail
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline
cd "$WD"
Rscript <script>.R
echo $(date +%T)
```

## Storage policy

Project roots:

- live: `/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline`
- ephemeral: `/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline`

Live storage must be sufficient to reproduce every final figure/table if
ephemeral storage is deleted. Save to live when a file:

- is read by another script;
- contains per-cell scores, states, gene lists, enrichment, normalized matrices,
  assignments, source data, or model summaries;
- is needed to regenerate a figure without rerunning the producing analysis.

Ephemeral storage is allowed only for exceptionally large, reproducible caches
or tool work directories. When uncertain, save to both. Create output
directories with `recursive=TRUE, showWarnings=FALSE`.

Long workflows use `intermediate/`, `tables/`, `figures/`, `logs/`, and
`reports/`. Where practical support `PDO_FORCE_REBUILD=1` and
`PDO_REPLOT_ONLY=1`; run logs record inputs, parameters, cache reuse, outputs,
times, and relevant package/session versions.

## Analysis governance

Every new or substantially updated analysis script starts with a registry block
modeled on the scRef centred scripts. It states:

- status and exact script path;
- methodology path, or `none` for a simple direct export/plot;
- `analysis/ANALYSIS_MAP.md`;
- accurate description;
- exact inputs and outputs;
- downstream use;
- cache/replot behavior;
- run command and conda environment.

Complex workflows require a mirrored methodology file detailed enough to
reconstruct the method and understand every consequential threshold, reference,
statistical unit, limitation, and storage decision. Simple deterministic
exports/plots do not need a standalone methodology. Retain obsolete methodology
only as clearly marked legacy provenance.

Update `analysis/ANALYSIS_MAP.md` whenever a script is added, renamed,
superseded, or gains a new dependency. Update this file only for durable
repository-wide rules, canonical objects, external dependencies, or recurring
technical constraints—not for per-script history or results.

New code sources shared constants/helpers from:

- `analysis/shared/Auto_pdo_analysis_config.R`
- `analysis/shared/Auto_pdo_analysis_helpers.R`

## Current centred MP/state route

Canonical live inputs:

| Object | Path under `PDOs_outs/` |
| :--- | :--- |
| post-QC list | `PDOs_list_PDOs.rds` |
| merged Seurat | `PDOs_merged.rds` |
| final MP genes | `centred_mp_refinement/merged_refined_mp_genes.rds` |
| final MP weights | `centred_mp_refinement/merged_refined_mp_gene_weights.rds` |
| final UCell | `centred_mp_refinement/merged_refined_ucell_scores.rds` |
| final states | `centred_mp_refinement/centred_refined_noreg_states.rds` |
| adjusted MP matrix | `centred_mp_refinement/centred_refined_noreg_mp_adj.rds` |
| group maxima | `centred_mp_refinement/centred_refined_noreg_group_max.rds` |
| MP grouping/order | `centred_mp_refinement/tables/centred_refined_mp_state_grouping.csv`; `centred_refined_mp_strict_order.rds` |
| enrichment | `centred_mp_refinement/cluster_enrich_centred.rds` |
| 3CA UCell | `UCell_3CA_MPs.rds` |

Do not use `MP_outs_default.rds`, `UCell_scores_filtered.rds`, uncentred
`geneNMF_metaprograms_nMP_13.rds`, `Auto_PDO_states_noreg.rds`, or
`Auto_PDO_final_states.rds` for new downstream work.

Centred QC:

- Parent MPs: keep at silhouette ≥0.2, split at 0<silhouette<0.2, remove at
  silhouette<0; keep/split also require coverage in at least 3 samples and more
  than 5 genes.
- Final refined MPs: coverage in at least 3 samples and at least 5 genes.
- Step 04 is the only final MP filter. Step 05 must not add a second threshold.
- No PDO-specific manual MP exclusion.
- State calls: unresolved if best normalized group score <0.5; hybrid if the
  best-minus-second score is <0.3; cell-cycle MPs never define a state.
- Named state vectors must preserve exact cell-barcode names through coercion.

## Canonical state and MP display order

States:

1. Classic proliferation — `#E41A1C`
2. Columnar-to-intestinal — `#4DAF4A`
3. Glandular differentiation — `#FF7F00`
4. Stress-adaptive — `#984EA3`
5. ECM-remodelling — `#A65628`
6. Motile-cilia differentiation — `#F781BF`

MP order follows those states, with cell-cycle MPs first:

| Group | MPs |
| :--- | :--- |
| Cell cycle/QC only | MP1, MP2, MP11, MP3 |
| Classic proliferation | MP19+ |
| Columnar-to-intestinal | MP14b, MP13b, MP5+, MP12, MP15 |
| Glandular differentiation | MP17+, MP8+ |
| Stress-adaptive | MP16b |
| ECM-remodelling | MP9 |
| Motile-cilia differentiation | MP18 |

Descriptions:

- MP1 G2/M cell cycle; MP2 G1/S cell cycle; MP11 Single-nucleus-associated cell cycle;
  MP3 Replication-dependent histones
- MP19+ MYC-associated proliferation
- MP14b Proliferative epithelial plasticity; MP13b Metabolic-detox columnar epithelium;
  MP5+ Inflammatory-reactive columnar epithelium; MP12 KRAS-active columnar epithelium;
  MP15 Intestinal metaplasia
- MP17+ Ciliated progenitor epithelium; MP8+ Secretory-transport glandular epithelium
- MP16b EMT/KRAS adaptive plasticity
- MP9 ECM-remodelling epithelium
- MP18 Motile-cilia differentiation

## Code conventions

R:

- `library()` calls at the top; no `require()`.
- snake_case variables/functions; preserve local style in old files.
- per-sample scripts read `commandArgs(trailingOnly=TRUE)`.
- use Seurat v5 `layer=` access with a compatibility fallback only where needed.
- use guard `stop()` calls for missing/empty critical inputs.
- save RDS/CSV/PDF/PNG with explicit, auditable paths.
- presentation figures must remain readable at exported slide dimensions.

Shell:

- `#!/bin/bash` first;
- timestamps at start/end;
- `module purge`, `module load tools/dev`, then conda initialization;
- quote paths/variables and use non-interactive commands.

## External references and compatibility

- 3CA MPs:
  `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`
- developmental references:
  `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/*.rds`
- cell-cycle genes:
  `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`
- oesophageal organoid differential-dependency reference:
  `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/CRISPR_screen_results_PDO.xlsx`
- scATLAS centred MP genes for cross-dataset comparisons:
  `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_genes.rds`
- scATLAS centred MP weights:
  `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/mp_refinement/intermediate/merged_refined_mp_gene_weights.rds`
- scATLAS epithelial object and current centred states:
  `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds`;
  `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/state_definition/intermediate/centred_refined_noreg_states.rds`
- scATLAS ranked centred state markers for cross-dataset comparisons:
  `/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/Metaprogrammes_Results/centred/state_markers/Auto_five_state_markers_ranked.csv`
- `dmtcp`: Seurat 5.1.0 with ggplot2 3.5.2; ggplot2 4.x breaks `VlnPlot()`.
- `gnmf`: Seurat 5.4.0/SeuratObject 5.3.0 with ggplot2 4.0.1 is compatible.
- On CX3, avoid `module load tools/bioinf` for Demuxafy/Souporcell.
  Singularity is available after `module purge`. Load `tools/prod` before
  BCFtools 1.22 or SAMtools 1.22.1.

Google Drive uploads always target `gdrive:IMPERIAL/`:

```bash
module load rclone
rclone copy <local_file> gdrive:IMPERIAL/ --progress
```

## Subagent model tier policy

If delegation is explicitly requested, restrict subagent models by the primary
model tier:

- Free primary: Free models only.
- 0.33X primary: Free or 0.33X models.
- Paid primary: any model.

Free: `opencode/big-pickle`, `opencode/minimax-m2.5-free`,
`opencode/trinity-large-preview-free`, `github-copilot/gpt-4.1`,
`github-copilot/gpt-4o`, `github-copilot/gpt-5-mini`.

0.33X: `github-copilot/gemini-3-flash-preview`,
`github-copilot/claude-haiku-4.5`,
`github-copilot/gpt-5.1-codex-mini`,
`github-copilot/grok-code-fast-1`.
