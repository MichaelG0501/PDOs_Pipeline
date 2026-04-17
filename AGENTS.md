# AGENTS.md — PDOs_Pipeline

Single-cell RNA-seq QC and analysis pipeline for OAC (Oesophageal Adenocarcinoma) Patient-Derived Organoids (PDOs). Runs on Imperial College HPC (PBS Pro scheduler). All computation is in **R** with **bash** PBS wrappers. Purely malignant cells (organoids) — no cell-type annotation or CNA inference needed.

## Repository Structure

```
*.R              — Core pipeline scripts (steps 1-3), executed via Rscript
*_master.sh      — PBS orchestrators that loop over samples and submit per-sample jobs
N_<Step>.sh      — PBS job scripts for each pipeline step
analysis/        — Downstream analysis R scripts (created as needed)
PDOs_outs/       — All pipeline outputs (gitignored)
temp/            — PBS stdout/stderr logs
```

## Pipeline Execution Order

| Step | Shell Script        | R Script           | Scope       | Conda Env |
|------|--------------------|--------------------|-------------|-----------|
| 1    | `1_QC_Pipeline.sh` | `QC_Pipeline.R`    | All samples | dmtcp     |
| 2    | `2_master.sh` → `2_NMF.sh` | `NMF.R`    | Per-sample  | dmtcp     |
| 3    | `3_geneNMF.sh`     | `geneNMF.R`        | All samples | gnmf      |

**Note:** Step 2 (per-sample NMF) uses `multiNMF()` internally in Step 3 via GeneNMF. The separate NMF.R step runs classical NMF per sample — skip if only using GeneNMF workflow.

### Excluded Sample
**SUR843T3_PDO** is excluded from all NMF and downstream analysis (removed in geneNMF.R and annotation scripts).

## Build / Run / Test Commands

There is no build system, linter, or test suite. All execution is via PBS `qsub`.

```bash
# Submit QC pipeline
qsub 1_QC_Pipeline.sh

# Submit per-sample NMF jobs (master scripts handle throttling)
qsub 2_master.sh

# Submit geneNMF metaprogram analysis
qsub 3_geneNMF.sh

# Submit per-sample NMF job manually
qsub -v sample="SUR1070_Treated_PDO" -N SUR1070_Treated_PDO 2_NMF.sh

# Run R interactively (for analysis scripts or debugging)
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
Rscript analysis/find_optimal_nmf.R

# For GeneNMF / UCell scripts, use the gnmf environment instead
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf
```

## HPC & File Safety Rules

These rules are **mandatory** for any agent operating in this repo:

1. **Working directory**: All outputs go to `PDOs_outs/`. Never write outside project paths.
2. **Conda init**: Always run `eval "$(~/miniforge3/bin/conda shell.bash hook)"` before activating envs.
3. **Interactive first**: Tasks under 8 cores / 64 GB → write only the `.R` script, no `.sh` wrapper. User runs interactively.
4. **PBS required**: Heavy tasks → must create PBS `.sh` script with `#PBS` resource headers.
5. **File naming**: New persistent files MUST be prefixed with `Auto_` (e.g., `Auto_analysis.R`).
6. **Modifying existing files**: New code MUST be wrapped in 20-hash comment blocks:
   ```r
   ####################
   # your new code here
   ####################
   ```
7. **No deleting/modifying** existing lines outside 20-hash blocks without permission.
8. **Test scripts**: Name `delete_<desc>.R` and delete immediately after use.
9. **Max concurrent PBS jobs**: 46 (throttled via `while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]`).

### PBS Job Template
```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
echo $(date +%T)
module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD
Rscript <script>.R
echo $(date +%T)
```

## Code Style Guidelines

### R Scripts

**Imports**: `library()` calls at top of file, one per line. Use double quotes for package names in `library("pkg")` or no quotes `library(pkg)` — both are used; match the file you're editing. No `require()`.

**Working directory**: Each script calls `setwd()` near the top to set the output directory (typically `PDOs_outs/`). All file paths are then relative to that.

**CLI arguments**: Per-sample scripts receive the sample name via `commandArgs(trailingOnly = TRUE)`:
```r
args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
```

**Variable naming**: `snake_case` for variables and functions. Short descriptive names. Examples:
- `pdos.list`, `pdos` — Seurat objects (PDO-specific)
- `filtered`, `merged_obj` — pipeline intermediates
- `geneNMF.programs`, `geneNMF.metaprograms` — GeneNMF outputs
- `mp.genes`, `mp_gene_lists` — metaprogram gene collections

**Data structures**: Heavy use of Seurat objects (`CreateSeuratObject`, metadata in `@meta.data`, assay layers `counts`, `data`, `CPM`). Lists of Seurat objects keyed by sample name.

**Tidyverse style**: Pipe-heavy with `%>%` and `|>`. `dplyr` verbs, `tidyr::pivot_longer/wider`, `purrr::map/imap`.

**Plotting**: `ggplot2` with `theme_minimal()` or `theme_classic()`. Plots saved via `ggsave()` or `pdf()`/`png()` + `dev.off()`. Composite layouts with `patchwork` (`|`, `/`) or `gridExtra::grid.arrange()`.

**Parallelism in R**: `mclapply()` from `parallel` package, or `foreach`/`doParallel`. NMF uses `nmf.options(parallel = 6)`.

**Output format**: `.rds` files for data, `.png` / `.pdf` for plots, `.csv` for summary tables. Saved to `PDOs_outs/` or `PDOs_outs/by_samples/<sample>/`.

### File Organization Pattern

Per-sample outputs follow: `PDOs_outs/by_samples/<SUR_ID>_<Condition>_PDO/`
- `<sample>.rds` — post-QC Seurat object
- `<sample>_rank4_9_nrun10.RDS` — NMF results
- Sentinel files: `no_cell` — skip markers

### Shell Scripts

**Shebang**: `#!/bin/bash` — always first line.

**PBS directives**: `#PBS -l select=...`, `#PBS -l walltime=...`, `#PBS -N <name>`. Resource sizing varies by step.

**Timestamps**: Every script prints `echo $(date +%T)` at start and end.

**Module loading**: `module purge` then `module load tools/dev` before conda.

**Master script pattern** (for per-sample parallelism):
```bash
for sample_folder in PDOs_outs/by_samples/*_PDO/; do
  while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]; do sleep 180; done
  sample=$(basename "$sample_folder")
  qsub -v sample=$sample -N $sample <step>.sh
done
```

### Error Handling

- Guard clauses with `stop()` for fatal conditions (e.g., not enough cells).
- Sentinel `.rds` files written before `stop()` to mark samples to skip in subsequent steps.
- `tryCatch()` for subsetting operations that might produce empty results.

### Key R Packages

**Core**: Seurat, dplyr, tidyr, purrr, ggplot2, readxl, Matrix, parallel
**Specialised**: NMF, GeneNMF, UCell, fgsea, msigdbr, DoubletFinder
**Plotting**: patchwork, gridExtra, RColorBrewer, circlize, scales
**Enrichment**: clusterProfiler, org.Hs.eg.db, enrichplot, pheatmap

## Key Shared Data Objects

All paths are relative to `PDOs_outs/` unless absolute paths are specified.

| Object | Path | Description |
| :--- | :--- | :--- |
| QC sample list | `PDOs_list_PDOs.rds` | Named list of post-QC Seurat objects per sample |
| Merged Seurat | `PDOs_merged.rds` | Merged Seurat object with PCA, UMAP, clustering |
| Metadata | `PDOs_all_meta.rds` | Per-cell metadata with clinical variables |
| Filtering summary | `filtering_summary_PDOs.csv` | Cell count per filtering step |
| GeneNMF programs | `geneNMF_outs.rds` | Raw multiNMF output object |
| MP results per nMP | `Metaprogrammes_Results/geneNMF_metaprograms_nMP_{k}.rds` | getMetaPrograms output for each nMP value |
| Optimal MP result | `MP_outs_default.rds` | Selected optimal nMP result |
| GO enrichment | `GO_outs.rds` | GSEA results per metaprogram |
| Final Seurat | `PDOs_final.rds` | Merged object with UCell MP scores |
| Enrichment results | `cluster_enrich.rds` | MP x database enrichment results |
| Surface marker ranked table | `Auto_five_state_surface_markers/Auto_five_state_surface_marker_ranked.csv` | Full PDO-only surface-marker scoring table with annotation support, topology flags, and target/off-target margins |
| Surface marker workbook | `Auto_five_state_surface_markers/Auto_five_state_surface_marker_candidates.xlsx` | Ranked FACS-oriented cell-surface candidate sheets across the five finalized PDO states |

## NMF Analysis Approach

### Finding Optimal Number of Metaprograms (nMP)

The pipeline uses the same approach as the scRef_Pipeline to determine optimal nMP:

1. **Run `multiNMF()`** once to get program-level NMF results (`geneNMF_outs.rds`)
2. **Loop `getMetaPrograms()`** over a range of nMP values (e.g., 4:20), saving each result as `geneNMF_metaprograms_nMP_{k}.rds`
3. **Evaluate** using:
   - **Silhouette analysis**: Average silhouette width from cosine distance matrix
   - **Elbow method (WSS)**: Within-cluster sum of squares based on cosine distance
4. **Select** the nMP at the inflection point (kneedle algorithm — maximum distance from diagonal)
5. **Filter** MPs with silhouette < 0 before downstream analysis

Key script: `analysis/metaprograms/Auto_find_optimal_nmf.R` — reads all `geneNMF_metaprograms_nMP_{k}.rds` files, applies kneedle inflection-point detection, and plots silhouette + WSS. **PDOs optimal: nMP=13** (silhouette inflection=13, WSS elbow=16).

## Enrichment Annotation

Metaprogram annotation uses enrichment analysis against multiple databases (same as scRef_Pipeline):
- **Hallmark** gene sets (MSigDB H)
- **GO Biological Process** (via clusterProfiler)
- **3CA metaprograms** (cancer cell atlas reference)
- **Developmental stages**: Early Embryogenesis, Organogenesis (major/sub), Normal Development (long/short), Adult Epithelium, Barretts Oesophagus

Reference data paths:
- 3CA MPs: `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`
- Developmental: `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/*.rds`

## PDO Nomenclature (States & Metaprograms)

Standardized nomenclature for PDO-specific analysis (Approach B):

| Cell State | Key Metaprograms | Color |
| :--- | :--- | :--- |
| **Classic Proliferative** | MP5 | `#E41A1C` (Red) |
| **Basal to Intest. Meta** | MP4 | `#4DAF4A` (Green) |
| **Stress-adaptive**       | MP10, MP9 | `#984EA3` (Purple) |
| **SMG-like Metaplasia**   | MP8 | `#FF7F00` (Orange) |
| **3CA_EMT_and_Protein_maturation** | 3CA_EMT_III, 3CA_ProtMat | `#377EB8` (Blue) |

**Metaprogram Descriptions:**

| MP | Standardized Description |
| :--- | :--- |
| **MP6** | G2M Cell Cycle |
| **MP7** | DNA repair |
| **MP5** | MYC-related Proliferation |
| **MP1** | G2M checkpoint |
| **MP3** | G1S Cell Cycle |
| **MP8** | Columnar Progenitor |
| **MP10** | Inflammatory Stress Epi. |
| **MP9** | ECM Remodeling Epi. |
| **MP4** | Intestinal Metaplasia |

Note: "Epi." is used instead of "(Epi.)" for consistency with scRef standards.

## Critical Recurring Patterns

**MP Silhouette Filtering**
Filter out MPs with a silhouette score below 0 before any downstream analysis. This is a strict requirement.
```r
bad_mps <- which(geneNMF.metaprograms$metaprograms.metrics$silhouette < 0)
bad_mp_names <- paste0("MP", bad_mps)
mp.genes <- mp.genes[!names(mp.genes) %in% bad_mp_names]
```

**MP Sample-Coverage Filtering (PDO-specific update)**
For PDO downstream annotation/correlation, also remove MPs with sample coverage < 25% (i.e., < 5 of 20 samples). This removes sparse MPs (currently MP11-13).
```r
coverage_tbl <- geneNMF.metaprograms$metaprograms.metrics$sampleCoverage
names(coverage_tbl) <- paste0("MP", seq_along(coverage_tbl))
low_coverage_mps <- names(coverage_tbl)[coverage_tbl < 0.25]
mp.genes <- mp.genes[!names(mp.genes) %in% low_coverage_mps]
```

**PDO Enrichment Display Order**
For PDO enrichment heatmaps, order metaprograms using reversed tree order:
```r
tree_order <- geneNMF.metaprograms$programs.tree$order
ordered_clusters <- geneNMF.metaprograms$programs.clusters[tree_order]
mp_tree_order <- rev(unique(ordered_clusters))
```

**Sample Exclusion**
SUR843T3_PDO is excluded from all NMF and downstream analysis:
```r
pdos.list$SUR843T3_PDO <- NULL
pdos <- subset(pdos, subset = orig.ident != "SUR843T3_PDO")
```

**PDO Surface-Marker Prioritization**
For FACS candidate ranking, reuse `Auto_five_state_marker_summary.csv` markers with `hit_sample_n > 0`, `best_state_match`, and `specificity_gap > 0`. Then require cell-surface annotation support, extracellular topology support, membrane-anchor or surfaceome support, and positive target-vs-off-target margins in both cell prevalence and mean expression.
```r
supported_for_facs <- annotation_surface_flag &
  extracellular_epitope_flag &
  (membrane_anchor_flag | surfaceome_surface_flag) &
  pct_margin > 0 &
  expr_margin > 0

recommended_for_facs <- supported_for_facs &
  (surfaceome_surface_flag | (uniprot_surface_flag & membrane_anchor_flag))
```

**Environment Usage**
- `gnmf` conda env: Use for UCell scoring and GeneNMF scripts.
- `dmtcp` conda env: Use for general Seurat and analysis tasks.

### Git Conventions

- `.gitignore`: ignores `PDOs_outs/` and PBS job output files (`*.[oe]*`)
- Commit messages: short, lowercase descriptions
- No tests, no CI/CD

## Subagent Model Tier Policy (MANDATORY)

When delegating work to subagents or background tasks, you **MUST** restrict model choices based on the tier of the primary model you are currently running as. Check your own model identity and apply the corresponding rule:

### Tier Classification

**Free Tier** (zero cost):
- `opencode/big-pickle`
- `opencode/minimax-m2.5-free`
- `opencode/trinity-large-preview-free`
- `github-copilot/gpt-4.1`
- `github-copilot/gpt-4o`
- `github-copilot/gpt-5-mini`

**0.33X Tier** (reduced cost):
- `github-copilot/gemini-3-flash-preview`
- `github-copilot/claude-haiku-4.5`
- `github-copilot/gpt-5.1-codex-mini`
- `github-copilot/grok-code-fast-1`

**All Other Models** = Paid Tier (full cost)

### Delegation Rules

| Your Primary Model Tier | Allowed Subagent Models |
|------------------------|------------------------|
| **Free** | Free tier models ONLY |
| **0.33X** | Free + 0.33X tier models |
| **Paid** (any other) | Any model — no restrictions |

## Uploading Files to Google Drive (rclone)

rclone is configured with a remote named `gdrive`. Always upload into the `IMPERIAL/` folder:

```bash
module load rclone
rclone copy <local_file> gdrive:IMPERIAL/ --progress
```

## AGENTS.md Living Document Rule

Future agents must update this file when they:
- Find a new analysis script and define its purpose.
- Locate new input or output file paths.
- Spot a recurring pattern or technical hurdle.
- Create a new `Auto_` script.

Append new findings to the appropriate section. Don't rewrite existing documentation unless fixing an error.

## Discovered Compatibility Notes

### ggplot2 + Seurat Compatibility
- **dmtcp env**: ggplot2 was downgraded from 4.0.2 → 3.5.2 because Seurat 5.1.0's `VlnPlot()` crashes internally with ggplot2 4.0.x (`aes_string()` deprecation + `+` operator changes). survminer 0.5.1 and ggpubr remain compatible with ggplot2 3.5.2.
- **gnmf env**: ggplot2 4.0.1 + Seurat 5.4.0 + SeuratObject 5.3.0 — **COMPATIBLE**. Seurat 5.4.0 works with ggplot2 4.0.x (updated internals). No downgrade needed.

### Analysis Directory Structure

```
analysis/
  cell_states/     — Auto_states_topmp_hybrid.R (approach B, PDO-adapted), Auto_PDO_state_concordance.R (scRef MPs vs PDO MPs concordance), Auto_pdo_overall_state_proportions.R (overall proportions barplot), Auto_sample_abundance_pdo.R (sample abundance with clinical annotations), Auto_marker_comparison_excel.R (cross-dataset scATLAS vs PDO marker comparison Excel), Auto_five_state_surface_markers.R (PDO-only FACS surface-marker prioritization from five-state marker outputs)
  clinical/        — Auto_clinical_mp_ucell_plots.R, Auto_clinical_variable_plots.R, Auto_survival_clinical_mps.R
  cnv/             — CNV_filter.R, cnv_profile.R, plot_CNV.R (copied from live PDOs)
  enrichment/      — Auto_enrichment_annotation.R, enrichment_extract.R, enrichment_plotting.R, enrich_plot.R, scGSEA.R, wnt_enrich.R
  methodology/     — Auto_five_state_surface_marker_methodology.md (operational description of the PDO FACS surface-marker workflow)
  metaprograms/    — Auto_find_optimal_nmf.R, Auto_extend_nMP_range.R, Auto_update_optimal_mp.R, Auto_ucell_vlnplot.R, Auto_mp_correlation_pdo.R, PDO_mp_correlation_crossdata.R, MP_analysis_pdos.R, robust_NMF.R, nmf_plot.R, Find_NMF.R, MP_dist.R, mp_ucell_scoring.R, robust_nmf_scref.R
  plotting/        — heatmap.R
```

### Auto_ Script Dependencies

| Script | Input Dependencies | Conda Env |
| :--- | :--- | :--- |
| `Auto_find_optimal_nmf.R` | `Metaprogrammes_Results/geneNMF_metaprograms_nMP_{4..35}.rds` | dmtcp |
| `Auto_extend_nMP_range.R` | `geneNMF_outs.rds` | gnmf |
| `Auto_update_optimal_mp.R` | `MP_outs_default.rds`, `PDOs_merged.rds` | gnmf |
| `Auto_ucell_vlnplot.R` | `MP_outs_default.rds`, `PDOs_merged.rds` | gnmf |
| `Auto_enrichment_annotation.R` | `MP_outs_default.rds` (optimal nMP result) | dmtcp |
| `Auto_mp_correlation_pdo.R` | `PDOs_final.rds`, `MP_outs_default.rds` | dmtcp |
| `Auto_states_topmp_hybrid.R` | `PDOs_merged.rds`, `MP_outs_default.rds`, `UCell_scores_filtered.rds` | dmtcp |
| `Auto_clinical_mp_ucell_plots.R` | `PDOs_merged.rds`, `MP_outs_default.rds`, `UCell_scores_filtered.rds` | dmtcp |
| `Auto_sample_abundance_pdo.R` | `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `UCell_scores_filtered.rds`, Clinical Excel | dmtcp |
| `Auto_clinical_variable_plots.R` | `PDOs_merged.rds` | dmtcp |
| `Auto_survival_clinical_mps.R` | `PDOs_merged.rds`, `UCell_scores_filtered.rds` | dmtcp |
| `PDO_mp_correlation_crossdata.R` | `geneNMF_metaprograms_nMP_13.rds`, `UCell_3CA_MPs.rds`, `scRef/...` | dmtcp |
| `Auto_PDO_state_concordance.R` | `PDOs_final.rds`, `scRef/.../geneNMF_metaprograms_nMP_19.rds`, `Auto_PDO_states_noreg.rds` | dmtcp |
| `Auto_PDO_final_mp_scenic.R` | `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `UCell_scores_filtered.rds`, `UCell_3CA_MPs.rds`, `geneNMF_metaprograms_nMP_13.rds`, `Auto_PDO_unresolved_relabel_mp_coverage.csv`, `New_NMFs.csv` | dmtcp |
| `Auto_marker_comparison_excel.R` | scRef `Auto_six_state_markers_ranked.csv` + `state_specificity.rds`, PDO `Auto_five_state_markers_ranked.csv` + `state_specificity.rds` | dmtcp |
| `Auto_five_state_surface_markers.R` | `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `Auto_five_state_marker_summary.csv`, `Auto_five_state_markers_ranked.csv`, UniProt reviewed-human surface/topology TSV (download/cache), ETH surfaceome Table S3 workbook (download/cache) | dmtcp |

### Additional Analysis Scripts

- `analysis/cell_states/Auto_pdo_sn_matched_pair_comparison.R` — compares the matched PDO/snRNA-seq pairs `SUR680T3_PDO -> H_post_T1_biopsy` and `SUR791T3_PDO -> L_post_T1_biopsy` using finalized state proportions plus modality-specific top-metaprogram proportions; writes a compact comparison report and CSV summaries to `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/`.
- `analysis/cell_states/Auto_pdo_flot_matched_response.R` — analyses the four matched untreated/FLOT-treated PDO pairs (`SUR1070`, `SUR1090`, `SUR1072`, `SUR1181`) using finalized state abundance shifts, paired pseudobulk Hallmark-response deltas, heuristic state-fate summaries, and paired edgeR state-specific pseudobulk DE; writes presentation-ready figures plus DEG tables to `PDOs_outs/Auto_pdo_flot_matched_response/`.
- `analysis/cell_states/Auto_marker_comparison_excel.R` — canonical cross-dataset marker-comparison workflow for scATLAS vs PDO; includes the combined 3-page marker heatmaps and the top-5 workbook output. The root `Auto_append.R` fragment is redundant.
- `analysis/cell_states/Auto_PDO_scAtlas_scenic_comparison.R` — canonical SCENIC comparison workflow for scATLAS vs PDO; includes RSS gap calculation, the 3-page RSS heatmaps, and the separate top-5 workbook output. The root `Auto_append_scenic.R`, `Auto_fix_rss_gap.R`, and `Auto_update_scenic.R` fragments are redundant.
- `analysis/cell_states/Auto_append_marker_excel.R` — legacy helper fragment for the top-5 workbook layout; its logic is already merged into `Auto_marker_comparison_excel.R`.

### Additional Auto_ Script Dependencies

- `Auto_pdo_sn_matched_pair_comparison.R`
  Inputs: PDO `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `UCell_scores_filtered.rds`, `Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`; snRNA-seq `/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs/snSeq_malignant_epi.rds`, `Auto_final_states.rds` (fallback `Auto_topmp_v2_noreg_states_B.rds`), `Metaprogrammes_Results/UCell_nMP19_filtered.rds`, `Metaprogrammes_Results/geneNMF_metaprograms_nMP_19.rds`
  Env: `dmtcp`
- `Auto_pdo_flot_matched_response.R`
  Inputs: `PDOs_merged.rds`, `Auto_PDO_final_states.rds`; Hallmark gene sets via `msigdbr`; matched sample IDs `SUR1070/SUR1090/SUR1072/SUR1181` treated vs untreated. Uses RNA `counts` for paired pseudobulk edgeR and RNA `data` for support-score UMAP overlays.
  Env: `dmtcp`
  Notes: paired pseudobulk analyses require at least 20 cells per sample-state; `3CA_EMT_and_Protein_maturation` retains only 3 valid matched pairs because `SUR1070_Treated_PDO` falls below threshold.

### Additional Output Paths

- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_comparison.pdf` — combined matched-pair figure with finalized-state bars and modality-specific top-MP proportion panels.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_comparison.png` — high-resolution PNG version of the matched-pair comparison figure.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_state_proportions.csv` — per-pair finalized state proportion table for the two PDO/snRNA-seq matches.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_mp_proportions.csv` — per-pair top-metaprogram proportion table for each matched PDO and snRNA-seq sample.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_summary.csv` — pair-level summary metrics including state overlap, correlations, dominant states, and dominant top MPs.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_summary_panels.pdf` — three-panel matched-FLOT summary figure combining paired state compositions, patient-level state log2 fold-changes, and paired pseudobulk composite response deltas.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_pathway_heatmap.pdf` — pathway-response heatmap with Hallmark treated-minus-untreated deltas for each patient-state combination.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_recurrent_deg_heatmaps.pdf` — one-page-per-state paired pseudobulk DEG heatmaps showing per-patient logFC for recurrent treated-vs-untreated genes.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_state_fate_summary.csv` — heuristic state-fate summary table combining abundance change, paired pseudobulk composite deltas, contributing pair counts, and interpretation calls.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_pseudobulk_pair_eligibility.csv` — per-state per-patient eligibility table for paired pseudobulk analyses under the 20-cell threshold.
- `PDOs_outs/Auto_pdo_flot_matched_response/pseudobulk_deg/Auto_pdo_flot_matched_deg_<state>.csv` — full paired edgeR state-specific pseudobulk DEG tables.

####################
### Additional Analysis Scripts

- `analysis/metaprograms/Auto_3CA_pseudobulk_correlation_crossdata.R` — compares pan-cancer 3CA metaprogram scores between scATLAS primary tumour and two PDO datasets side-by-side: OAC PDO single-cell pseudobulk vs scATLAS pseudobulk, and OSCC/ESCC PDO bulk RNA-seq (GSE269447 tumor organoids only) vs scATLAS pseudobulk. Scores all datasets with shared 3CA gene sets using `ScoreSignatures_UCell()` on bulk-like sample matrices and labels the paired scatter plots with Spearman statistics.

### Additional Auto_ Script Dependencies

- `Auto_3CA_pseudobulk_correlation_crossdata.R`
  Inputs: `PDOs_merged.rds`; scATLAS `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/EAC_Ref_epi.rds`; 3CA gene sets `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv`; OSCC bulk download/extract staging `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/raw_txt/GSM*_Tumor-Org_TPM.txt.gz`
  Env: `dmtcp`
  Notes: excludes `SUR843T3_PDO`; uses sample-level pseudobulk counts for OAC PDO and scATLAS, uses the raw count column from each `GSE269447` tumor-organoid quantification file for OSCC, collapses Ensembl IDs to gene symbols with `org.Hs.eg.db`, and keeps only 3CA MPs with at least 5 genes present in every dataset before scoring.

### Additional Output Paths

- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata.pdf` — side-by-side scatter plot for `OAC PDO vs scATLAS` and `OSCC PDO vs scATLAS` with shared-axis pseudobulk/bulk 3CA UCell scores and Spearman annotations.
- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata.png` — PNG version of the side-by-side pseudobulk/bulk 3CA correlation plot.
- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata_summary.csv` — per-comparison summary table with Spearman rho, p-value, number of shared MPs, highlighted MPs, and target/scATLAS sample counts.
- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata_mean_scores.csv` — dataset-level mean 3CA scores for each shared MP across scATLAS, OAC PDO, and OSCC PDO.
- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_correlation_crossdata_gene_overlap.csv` — per-dataset gene-overlap table reporting the number and fraction of genes retained for each shared 3CA MP after symbol harmonization.
- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_scores_pdo.csv` — OAC PDO sample-by-MP pseudobulk UCell score matrix.
- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_scores_scatlas.csv` — scATLAS sample-by-MP pseudobulk UCell score matrix.
- `PDOs_outs/Auto_3CA_pseudobulk_correlation_crossdata/Auto_3CA_pseudobulk_scores_oscc.csv` — OSCC PDO bulk sample-by-MP UCell score matrix.

### Additional External Data Paths

- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_OSCC_PDO_GSE269447/` — external staging directory for GEO `GSE269447` OSCC/ESCC PDO bulk RNA-seq download. Contains `GSE269447_RAW.tar` plus extracted `raw_txt/GSM*_Tumor-Org_TPM.txt.gz` files used by `Auto_3CA_pseudobulk_correlation_crossdata.R`.
####################
