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
5. **Live Logging**: Always use live streaming log file mode by adding `#PBS -koed` to the submission script. This ensures standard out and standard error are written to their final destination as the job is running, allowing for real-time monitoring from login nodes.
6. **File naming**: New persistent files MUST be prefixed with `Auto_` (e.g., `Auto_analysis.R`).
7. **Modifying existing files**: New code MUST be wrapped in 20-hash comment blocks:
   ```r
   ####################
   # your new code here
   ####################
   ```
8. **No deleting/modifying** existing lines outside 20-hash blocks without permission.
9. **Test scripts**: Name `delete_<desc>.R` and delete immediately after use.
10. **Max concurrent PBS jobs**: 46 (throttled via `while [[ $(qstat | grep sg3723 | wc -l) -gt 46 ]]`).

### PBS Job Template
```bash
#!/bin/bash
#PBS -l select=1:ncpus=<N>:mem=<M>gb
#PBS -l walltime=<HH:MM:SS>
#PBS -N <jobname>
#PBS -koed
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
| MP UCell scores | `UCell_scores_filtered.rds` | Metaprogram scores for the filtered MP set |
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
| **SMG-like Metaplasia**   | MP8 | `#FF7F00` (Orange) |
| **Stress-adaptive**       | MP10, MP9 | `#984EA3` (Purple) |
| **3CA_EMT_and_Protein_maturation** | 3CA_EMT_III, 3CA_ProtMat | `#377EB8` (Blue) |

### Canonical State And MP Display Order

For all new code and any future updates touching state order, display finalized malignant states in this order: `Classic Proliferative`, `Basal to Intest. Meta`, `SMG-like Metaplasia`, `Stress-adaptive`, optional `Immune Infiltrating` for scATLAS-style state sets when present, then `3CA_EMT_and_Protein_maturation`.

Metaprogram display order must follow this state order. Within each state, order MPs by the current `mp_tree_order` rather than by numeric MP ID when `mp_tree_order` is available.

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

**Final-State Vector Normalization**
When normalizing or remapping a saved per-cell state vector such as `Auto_PDO_final_states.rds`, preserve the original cell-barcode names. Coercing the vector to plain character without restoring `names()` will silently break downstream lookups like `state_vec[Cells(seurat_obj)]`.
```r
state_names <- names(state_vec)
state_vec <- as.character(state_vec)
names(state_vec) <- state_names
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
  metaprograms/    — Auto_find_optimal_nmf.R, Auto_extend_nMP_range.R, Auto_update_optimal_mp.R, Auto_mp_correlation_pdo.R, PDO_mp_correlation_crossdata.R, MP_analysis_pdos.R, robust_NMF.R, nmf_plot.R, Find_NMF.R, MP_dist.R, mp_ucell_scoring.R, robust_nmf_scref.R
  plotting/        — heatmap.R
```

### Auto_ Script Dependencies

| Script | Input Dependencies | Conda Env |
| :--- | :--- | :--- |
| `Auto_find_optimal_nmf.R` | `Metaprogrammes_Results/geneNMF_metaprograms_nMP_{4..35}.rds` | dmtcp |
| `Auto_extend_nMP_range.R` | `geneNMF_outs.rds` | gnmf |
| `Auto_update_optimal_mp.R` | `MP_outs_default.rds`, `PDOs_merged.rds` | gnmf |
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
- `analysis/cell_states/Auto_pdo_sn_matched_pair_comparison.R` (update) — now writes a two-page PDF: page 1 focuses on non-cell-cycle state-defining top MPs and renormalized stacked state bars excluding `Unresolved`/`Hybrid`; page 2 keeps the full-view comparison but without subtitles, with centered legends to avoid left-edge clipping.
- `analysis/cell_states/Auto_pdo_flot_matched_response.R` — analyses the four matched untreated/FLOT-treated PDO pairs (`SUR1070`, `SUR1090`, `SUR1072`, `SUR1181`) using finalized state abundance shifts, paired pseudobulk Hallmark-response deltas, heuristic state-fate summaries, and paired edgeR state-specific pseudobulk DE; writes presentation-ready figures plus DEG tables to `PDOs_outs/Auto_pdo_flot_matched_response/`.
- `analysis/cell_states/Auto_pdo_flot_presentation_final.R` — consumes the cached matched-FLOT response outputs and rebuilds the final multi-page presentation PDF, including the pathway matrix, composite response, state abundance, MP change, hybrid-only abundance, recurrent DEG, fate summary, and UMAP support pages.
- `analysis/cell_states/PDO_finalize_states.R` — regenerates `Auto_PDO_final_states.rds` from `unresolved_states/Auto_PDO_unresolved_relabel_states.rds`, merges the stray 3CA respiration/cell-cycle and EMT/protein-maturation labels into the canonical finalized PDO states, and redraws the per-cell heatmap plus state-proportion summaries before the downstream TCGA survival-volcano step.
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
- `Auto_pdo_flot_presentation_final.R`
  Inputs: `Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_response_results.rds`, `PDOs_merged.rds`, `Auto_PDO_final_states.rds`; Hallmark gene sets via `msigdbr` for the support-score UMAP page.
  Env: `dmtcp`
- `PDO_finalize_states.R`
  Inputs: `PDOs_merged.rds`, `UCell_scores_filtered.rds`, `UCell_3CA_MPs.rds`, `Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`, `unresolved_states/Auto_PDO_unresolved_relabel_states.rds`; downstream TCGA volcano step also expects `/rds/general/project/spatialtranscriptomics/ephemeral/TCGA/INPUT/TCGA_ESCA_TPM_CIBERSORTx_Mixture.txt`.
  Env: `dmtcp`
  Notes: `Auto_PDO_final_states.rds` is written before the TCGA input is read, so rerunning this script still refreshes the finalized state vector even if the TCGA file is absent or unreadable.

### Additional Output Paths

- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_comparison.pdf` — combined matched-pair figure with finalized-state bars and modality-specific top-MP proportion panels.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_comparison.png` — high-resolution PNG version of the matched-pair comparison figure.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_state_proportions.csv` — per-pair finalized state proportion table for the two PDO/snRNA-seq matches.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_mp_proportions.csv` — per-pair top-metaprogram proportion table for each matched PDO and snRNA-seq sample.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_summary.csv` — pair-level summary metrics including state overlap, correlations, dominant states, and dominant top MPs.
- `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/Auto_pdo_sn_matched_pair_comparison.pdf` (update) — single two-page PDF export only (`.png` export removed): page 1 is the state-defining/non-cell-cycle focus view and page 2 is the full-view comparison without subtitles.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_summary_panels.pdf` — three-panel matched-FLOT summary figure combining paired state compositions, patient-level state log2 fold-changes, and paired pseudobulk composite response deltas.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_pathway_heatmap.pdf` — pathway-response heatmap with Hallmark treated-minus-untreated deltas for each patient-state combination.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_recurrent_deg_heatmaps.pdf` — one-page-per-state paired pseudobulk DEG heatmaps showing per-patient logFC for recurrent treated-vs-untreated genes.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_response_results.rds` — cached matched-FLOT analysis object consumed by `Auto_pdo_flot_presentation_final.R`.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_presentation_final.pdf` — final eight-page presentation export for the matched-FLOT analysis.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_state_fate_summary.csv` — heuristic state-fate summary table combining abundance change, paired pseudobulk composite deltas, contributing pair counts, and interpretation calls.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_pseudobulk_pair_eligibility.csv` — per-state per-patient eligibility table for paired pseudobulk analyses under the 20-cell threshold.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_mp_expression_changes.csv` — per-patient treated-vs-untreated log2FC summary for the finalized top metaprograms shown in the MP change panel.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_matched_hybrid_abundance_changes.csv` — per-patient treated-vs-untreated abundance shifts for the six pairwise hybrid classes built from the four non-3CA finalized states.
- `PDOs_outs/Auto_pdo_flot_matched_response/pseudobulk_deg/Auto_pdo_flot_matched_deg_<state>.csv` — full paired edgeR state-specific pseudobulk DEG tables.
- `PDOs_outs/Auto_PDO_final_states.rds` — finalized per-cell PDO state vector used by downstream state-abundance, FLOT-response, SCENIC, and marker workflows.

####################

####################
### Additional Analysis Scripts

- `analysis/cell_states/Auto_marker_selection_simulation.R` — validates selected PDO/scATLAS marker panels using two in silico experiments: marker-only single-cell assignment from RNA `data` layer marker gates based on marker-comparison target/off-target medians, and qRT-PCR-like paired shift simulations where one finalized PDO state increases while another decreases. Keeps immune-regulated markers as a PDO-negative-control panel but excludes that panel from assignable Simulation 1 states.

### Additional Auto_ Script Dependencies

- `Auto_marker_selection_simulation.R`
  Inputs: `Auto_five_state_markers/cache/pdos_state5_embedded.rds`, `Auto_five_state_markers/cache/state_specificity.rds` (fallback: `PDOs_merged.rds`, `Auto_PDO_final_states.rds`)
  Env: `dmtcp`
  Notes: excludes `SUR843T3_PDO`; samples only finalized five-state cells, so `Unresolved` in Simulation 1 is a below-gate marker call rather than an input state. Runtime/gates can be adjusted with `AUTO_MARKER_GATE_REPS`, `AUTO_MARKER_GATE_N`, `AUTO_MARKER_GATE_THRESHOLD_FRACTION`, `AUTO_MARKER_GATE_MIN_POSITIVE`, `AUTO_MARKER_QPCR_REPS`, `AUTO_MARKER_QPCR_N`, `AUTO_MARKER_QPCR_SHIFT_MIN`, and `AUTO_MARKER_QPCR_SHIFT_MAX`.

### Additional Output Paths

- `PDOs_outs/Auto_marker_selection_simulation/Auto_marker_panel_manifest.csv` — selected marker-panel manifest with PDO gene-availability flags.
- `analysis/methodology/Auto_marker_selection_simulation_methodology.md` — persistent methodology note for the marker-selection and qRT-PCR shift simulations.
- `PDOs_outs/Auto_marker_selection_simulation/Auto_marker_expression_validation_vs_marker_comparison.csv` — confirms selected marker target-state expression against the marker-comparison specificity cache.
- `PDOs_outs/Auto_marker_selection_simulation/Auto_marker_detection_by_state.csv` and `Auto_marker_detection_by_state_dotplot.pdf` — per-marker detection and sample-aware RNA `data` expression summaries by finalized PDO state, faceted by marker panel with colored state labels.
- `PDOs_outs/Auto_marker_selection_simulation/Auto_marker_gate_gene_positive_by_state.csv`, `Auto_marker_gate_panel_positive_by_state.csv`, `Auto_marker_gate_simulation_confusion_replicates.csv`, `Auto_marker_gate_simulation_metrics_replicates.csv`, `Auto_marker_gate_simulation_summary.csv`, and `Auto_marker_gate_simulation_confusion_heatmap.pdf` — marker-only single-cell gate diagnostics and assignment simulation outputs.
- `PDOs_outs/Auto_marker_selection_simulation/Auto_qpcr_abundance_simulation_replicates.csv`, `Auto_qpcr_gene_log2fc_replicates.csv`, `Auto_qpcr_abundance_simulation_summary.csv`, `Auto_qpcr_gene_log2fc_summary.csv`, `Auto_qpcr_shift_role_summary.csv`, `Auto_qpcr_transition_eval_replicates.csv`, `Auto_qpcr_transition_summary.csv`, `Auto_qpcr_example_shift_lineplots.pdf`, `Auto_qpcr_repeated_shift_line_summary.csv`, and `Auto_qpcr_repeated_shift_summary_lineplots.pdf` — qRT-PCR-like paired state-shift simulation outputs with condition-connected line plots.
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

####################
### Additional Analysis Scripts

- `analysis/cell_states/Auto_scATLAS_four_marker_specificity.R` — focused scATLAS-only marker-specificity visualization for four PDO-prioritized surface markers (`MUC13`, `CEACAM5`, `ROR1`, `PTPRG`) using the scATLAS specificity intermediate used by `Auto_marker_comparison_excel.R`; highlights target states (`Basal to Intestinal Metaplasia` for `MUC13/CEACAM5`, `SMG-like Metaplasia` for `ROR1/PTPRG`) and exports both figure and compact specificity summary table.

### Additional Auto_ Script Dependencies

- `Auto_scATLAS_four_marker_specificity.R`
  Inputs: scATLAS `ref_outs/Auto_six_state_markers/cache/state_specificity.rds` (from `/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/`)
  Env: `dmtcp`

### Additional Output Paths

- `PDOs_outs/Auto_four_marker_scATLAS_specificity/Auto_four_marker_scATLAS_specificity_heatmap.pdf` — publication-style four-marker x six-state scATLAS heatmap with target-state highlight circles.
- `PDOs_outs/Auto_four_marker_scATLAS_specificity/Auto_four_marker_scATLAS_specificity_heatmap.png` — high-resolution PNG version of the focused scATLAS specificity heatmap.
- `PDOs_outs/Auto_four_marker_scATLAS_specificity/Auto_four_marker_scATLAS_specificity_summary.csv` — per-marker target-state specificity summary including target expression, max off-target expression, margin, best-state match flag, and target rank.
####################

####################
### Additional Analysis Scripts

- `analysis/cell_states/Auto_drug_reversal_inputs.R` — prepares five-state PDO state-vs-rest drug-reversal inputs: pooled Seurat DEGs, ASGARD gene-list objects, scDrugPrio DEG tables, CLUE/CMap top-150 up/down GMT files, sparse matrix/metadata exports, and ASGARD anchor-gene diagnostics.
- `analysis/cell_states/Auto_drug_reversal_asgard.R` — ASGARD wrapper for negative-connectivity mono-drug reversal from the prepared PDO DEG lists; requires an ASGARD L1000 reference RDS or rankMatrix/gene/drug info paths supplied through `AUTO_ASGARD_*` environment variables.
- `analysis/cell_states/Auto_drug_reversal_scdrugprio.R` — scDrugPrio transcriptomic/network-prioritization wrapper using PDO DEGs, PPI resources, drug-target mappings, and optional pharmacological action annotations; deliberately does not run scDrug's cell-line cytotoxicity module.
- `analysis/cell_states/Auto_drug_reversal_clue_fallback.R` — direct CLUE/CMap L1000 Touchstone fallback that submits the top-150 up/down PDO signatures when `CLUE_API_KEY` or `CLUE_KEY` is available; otherwise leaves submission-ready GMT files and a status file.
- `analysis/cell_states/Auto_drug_reversal_consensus_visuals.R` — intersects top-100 ASGARD and top-100 scDrugPrio/fallback candidates per state, then exports overlap summaries, Venn diagrams, rank-rank scatter plots, and top-consensus drug-target expression dot plots.
- `analysis/cell_states/Auto_drug_reversal_method_visuals.R` — generates presentation-level method evidence plots for ASGARD, scDrugPrio, and local CMap/L1000: top-rank heatmaps, three-method top-100 overlap bars, final-overlap rank-support matrix, repeated target/MoA summaries, and local L1000 state-up/state-down signature profile plots.
- `analysis/cell_states/Auto_drug_reversal_predicted_reversion_visuals.R` — generates prediction-specific inhibitor evidence plots: per-method rank waterfall plots, LINCS/ASGARD anti-correlation scatter plots comparing PDO state-vs-rest logFC to inferred opposing drug coordinates, predicted state-flipping heatmaps, and DrugBank/PPI targeted-hub overlays for annotated final-overlap drugs.
- `analysis/cell_states/Auto_drug_reversal_scdrugprio_visuals.R` — focused scDrugPrio visualization script with waterfall, L1000 anti-correlation, and refined PPI hub overlay. The refined network writes edge/node CSVs with target logFC, action type, and `target_direction` classifications (`Counteracts`, `Mimics`, `No DEG information`, `Non-directional/unknown`) so activator/inhibitor direction is visually auditable.
- `analysis/cell_states/Auto_prepare_asgard_reference.R` — builds ASGARD tissue-specific L1000 rank-matrix references from downloaded `GSE70138`/`GSE92742` files using `Asgard::PrepareReference()` and writes reusable ASGARD path exports for the drug-reversal wrapper.

### Additional Auto_ Script Dependencies

- `Auto_drug_reversal_inputs.R`
  Inputs: `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, optional `Auto_five_state_markers/cache/pdos_state5_embedded.rds`
  Env: `dmtcp`
  Notes: excludes `SUR843T3_PDO`, `Hybrid`, `Unresolved`, and non-finalized extra 3CA labels; default DE mode is pooled Seurat `FindMarkers()` state-vs-rest with `AUTO_DRUG_DEG_MIN_PCT=0.01`, top signatures controlled by `AUTO_DRUG_SIGNATURE_TOP_N` (default 150), and matrix export controlled by `AUTO_EXPORT_DRUG_MATRIX`. Fresh reruns should use `AUTO_FORCE_DRUG_DEGS=1` and `AUTO_DRUG_DEG_MODE=findmarkers`. The script now writes per-state DEG checkpoints under `PDOs_outs/Auto_drug_reversal/deg_checkpoints/`, falls back cleanly if `Auto_drug_reversal/cache/Auto_drug_reversal_state5.rds` is corrupted, atomically rewrites that cache, and slims the large marker-cache fallback with `DietSeurat()` before saving.
- `Auto_drug_reversal_asgard.R`
  Inputs: `Auto_drug_reversal/asgard_inputs/Auto_asgard_gene_list.rds`; plus either `AUTO_ASGARD_DRUG_REF_RDS` or `AUTO_ASGARD_DRUG_RESPONSE`, `AUTO_ASGARD_GENE_INFO`, and `AUTO_ASGARD_DRUG_INFO`
  Env: dedicated `PDOs_outs/Auto_drug_reversal/conda/Auto_drug_reversal` env from `Auto_setup_drug_reversal_env.sh`
- `Auto_prepare_asgard_reference.R`
  Inputs: uncompressed GEO L1000 files under `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000/plain/`
  Env: dedicated drug-reversal env preferred; default target tissue is `stomach` because ASGARD/LINCS metadata does not include an oesophagus primary site. Override with `AUTO_ASGARD_TISSUE`.
- `Auto_drug_reversal_scdrugprio.R`
  Inputs: `Auto_drug_reversal/scdrugprio_inputs/Auto_scdrugprio_deg_<state>.txt`, `AUTO_SCDRUGPRIO_PPI`, `AUTO_SCDRUGPRIO_DRUG_TARGETS`, optional `AUTO_SCDRUGPRIO_PHARMA_EFFECT`
  Env: dedicated `PDOs_outs/Auto_drug_reversal/conda/Auto_drug_reversal` env from `Auto_setup_drug_reversal_env.sh`
  Notes: runtime can be adjusted with `AUTO_SCDRUGPRIO_CORES`, `AUTO_SCDRUGPRIO_RANDOM_ITERATIONS`, `AUTO_SCDRUGPRIO_PADJ`, `AUTO_SCDRUGPRIO_MIN_ABS_LOGFC`, and `AUTO_SCDRUGPRIO_MAX_DISEASE_GENES`. scDrugPrio ranking must prioritize pharmacological direction before network proximity: `direction_call == "mimicking"` is excluded by default with `AUTO_SCDRUGPRIO_EXCLUDE_MIMICS=1`, and final selected hits require at least one counteracting DEG target with `AUTO_SCDRUGPRIO_REQUIRE_COUNTERACTION=1`. Drugs whose targets have no DEG/direction information should not be treated as significant reversal hits.
- `Auto_drug_reversal_clue_fallback.R`
  Inputs: `Auto_drug_reversal/clue_inputs/Auto_clue_up_entrez.gmt`, `Auto_drug_reversal/clue_inputs/Auto_clue_down_entrez.gmt`, and `CLUE_API_KEY` or `CLUE_KEY`
  Env: dedicated drug-reversal env preferred; falls back to `dmtcp` if the env is absent.
  Notes: saved CLUE response JSON is redacted because CLUE can echo the API key; batch five-state top-50 submissions returned CLUE server OOM, so use `AUTO_CLUE_KEEP_STATES`, `AUTO_CLUE_RUN_LABEL`, `AUTO_CLUE_TOOL_ID=sig_gutc_tool`, and `AUTO_CLUE_JOB_ID` for per-state submission/polling.
- `Auto_drug_reversal_consensus_visuals.R`
  Inputs: `Auto_drug_reversal/asgard/Auto_asgard_ranked_drugs.csv`, `Auto_drug_reversal/scdrugprio/Auto_scdrugprio_ranked_drugs.csv`; optional fallback `Auto_drug_reversal/clue_fallback/Auto_clue_ranked_drugs.csv`; `PDOs_merged.rds` and `Auto_PDO_final_states.rds` for target expression if cached state object is absent.
  Env: `dmtcp`
- `Auto_drug_reversal_method_visuals.R`
  Inputs: `Auto_drug_reversal/asgard/Auto_asgard_ranked_drugs.csv`, `Auto_drug_reversal/scdrugprio/Auto_scdrugprio_ranked_drugs.csv`, `Auto_drug_reversal/clue_fallback/Auto_clue_ranked_drugs.csv`, `Auto_drug_reversal/Auto_drug_reversal_signature_top150.csv`, and `Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.csv`
  Env: `dmtcp`
  Notes: uses the canonical state order with SMG-like Metaplasia before Stress-adaptive; the local L1000 profile is a rank-space diagnostic showing whether selected drugs separate target-state upregulated and downregulated genes in the perturbation reference.
- `Auto_drug_reversal_predicted_reversion_visuals.R`
  Inputs: `Auto_drug_reversal/Auto_drug_reversal_degs_all_states.csv.gz`, `Auto_drug_reversal/Auto_drug_reversal_signature_top150.csv`, all three ranked drug tables, ASGARD reference paths, scDrugPrio PPI, and DrugBank target/action table.
  Env: `dmtcp`
  Notes: all outputs are predicted/computational evidence, not treated PDO measurements. The heatmap column named `Predicted treatment (opposing LINCS coordinate)` is derived from the inverse centered LINCS rank coordinate and should be described as a reference-signature hypothesis.
- `Auto_drug_reversal_scdrugprio_visuals.R`
  Inputs: `Auto_drug_reversal/scdrugprio/Auto_scdrugprio_ranked_drugs.csv`, `Auto_drug_reversal/Auto_drug_reversal_degs_all_states.csv.gz`, `Auto_drug_reversal/Auto_drug_reversal_signature_top150.csv`, ASGARD reference paths, scDrugPrio PPI, and DrugBank target/action table.
  Env: `dmtcp`
  Notes: refined PPI plots mark drug nodes separately and keep the activator/inhibitor/unknown action labels visible. The plot should use `Auto_scdrugprio_direction_audit_all_drugs.csv` when available so all top network-proximity candidates can be inspected, while `Auto_scdrugprio_ranked_drugs.csv` remains the strict final table after scDrugPrio pharmacological-action filtering. A target with `No DEG information` should not be interpreted as either counteracting or mimicking without additional evidence, and should not drive final scDrugPrio hit selection.

### Additional Output Paths

- `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_degs_all_states.csv.gz` — pooled five-state state-vs-rest DEG table for drug reversal signatures.
- `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_signature_top150.csv` — top 150 upregulated and top 150 downregulated genes per finalized PDO state.
- `PDOs_outs/Auto_drug_reversal/asgard_inputs/Auto_asgard_gene_list.rds` and `Auto_asgard_deg_<state>.txt` — ASGARD-ready ranked DEG inputs.
- `PDOs_outs/Auto_drug_reversal/scdrugprio_inputs/Auto_scdrugprio_deg_<state>.txt` — scDrugPrio-ready state DEG inputs.
- `PDOs_outs/Auto_drug_reversal/clue_inputs/Auto_clue_up_entrez.gmt`, `Auto_clue_down_entrez.gmt`, `Auto_clue_up_symbols.gmt`, and `Auto_clue_down_symbols.gmt` — CLUE/CMap direct query gene sets.
- `PDOs_outs/Auto_drug_reversal/matrix/Auto_drug_reversal_counts.mtx`, `Auto_drug_reversal_features.tsv`, `Auto_drug_reversal_barcodes.tsv`, and `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_metadata.csv` — sparse matrix and metadata exports for external wrappers.
- `PDOs_outs/Auto_drug_reversal/asgard/Auto_asgard_ranked_drugs.csv` — standardized ASGARD state-drug ranking when ASGARD references are available.
- `PDOs_outs/Auto_drug_reversal/scdrugprio/Auto_scdrugprio_ranked_drugs.csv` — standardized scDrugPrio state-drug ranking when PPI/drug-target resources are available.
- `PDOs_outs/Auto_drug_reversal/clue_fallback/Auto_clue_ranked_drugs.csv` — standardized direct CLUE fallback ranking when a CLUE API key is available.
- `PDOs_outs/Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.csv` and `.sh` — generated ASGARD reference paths after `Auto_prepare_asgard_reference.R` completes.
- `PDOs_outs/Auto_drug_reversal/consensus/Auto_drug_reversal_top100_overlap_summary.csv` and `Auto_drug_reversal_consensus_drugs.csv` — per-state top-100 overlap counts and prioritized consensus drug table.
- `PDOs_outs/Auto_drug_reversal/consensus/Auto_drug_reversal_venn_top100.pdf`, `Auto_drug_reversal_rank_rank_scatter.pdf/.png`, and `Auto_drug_reversal_mechanism_target_dotplot.pdf/.png` — presentation-ready consensus visualizations.
- `PDOs_outs/Auto_drug_reversal/method_visuals/Auto_drug_reversal_method_rank_heatmap.pdf/.png`, `Auto_drug_reversal_three_method_overlap_barplot.pdf/.png`, `Auto_drug_reversal_final_overlap_rank_matrix.pdf/.png`, `Auto_drug_reversal_final_overlap_target_frequency.pdf/.png`, `Auto_drug_reversal_final_overlap_moa_barplot.pdf/.png`, and `Auto_drug_reversal_l1000_signature_reversal_profiles.pdf/.png` — method-level and final-overlap presentation figures explaining why inhibitors were selected.
- `PDOs_outs/Auto_drug_reversal/predicted_reversion_visuals/Auto_drug_reversal_reversion_waterfall_by_method.pdf/.png`, `Auto_drug_reversal_predicted_anticorrelation_scatter.pdf/.png`, `Auto_drug_reversal_predicted_state_flipping_heatmap.pdf/.png`, and `Auto_drug_reversal_ppi_targeted_hub_overlay.pdf/.png` — prediction-specific visual evidence that selected inhibitors oppose PDO state signatures and/or target annotated PPI hubs.
- `PDOs_outs/Auto_drug_reversal/scdrugprio_visuals/Auto_scdrugprio_ppi_hub_refined.pdf/.png`, `Auto_scdrugprio_ppi_hub_refined_nodes.csv`, `Auto_scdrugprio_ppi_hub_refined_edges.csv`, `Auto_scdrugprio_waterfall.pdf/.png`, and `Auto_scdrugprio_reversion_scatter.pdf` — focused scDrugPrio visual outputs with explicit target-action direction and target DEG status.
- `analysis/methodology/Auto_drug_reversal_methodology.md` — operational methodology for the ASGARD/scDrugPrio/CLUE consensus reversal workflow.

### Additional External Data Paths

- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000/` — ASGARD L1000 staging directory for downloaded GEO `GSE70138`/`GSE92742` raw `.gz` files, uncompressed `.txt`/`.gctx` files, and generated `DrugReference/` tissue rank matrices. Download job `2529425.pbs-7` and reference build job `2529426.pbs-7` were submitted on 2026-04-23; build depends on the package-env job `2529424.pbs-7` and download job completing successfully.
- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/ppi.txt` — scDrugPrio explicit PPI network resource; two-column Entrez-ID PPIN (`Protein_A`, `Protein_B`) from the scDrugPrio/Figshare input bundle.
- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/all_drug_targets_drug_bank.txt` — scDrugPrio/DrugBank explicit drug-target and pharmacological action table. In this local file, `target_organism` stores the action label (e.g. inhibitor/agonist/antagonist) and `drug_action` stores the organism; `Auto_drug_reversal_scdrugprio.R` and `Auto_drug_reversal_method_visuals.R` account for this header swap.
####################

####################
### Additional Analysis Scripts

- `analysis/cell_states/Auto_drug_reversal_local_cmap.R` — direct local CMap fallback that scores each five-state top-150 up/down PDO signature against the ASGARD tissue-specific `rankMatrix.txt`, then annotates overlapping compounds with CLUE `perts` endpoint `target` and `moa` metadata when `CLUE_API_KEY` or `CLUE_KEY` is available.

### Additional Auto_ Script Dependencies

- `Auto_drug_reversal_local_cmap.R`
  Inputs: `Auto_drug_reversal/Auto_drug_reversal_signature_top150.csv`, `Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.csv`, and the referenced `stomach_rankMatrix.txt`, `stomach_gene_info.txt`, `stomach_drug_info.txt`
  Env: `dmtcp`
  Notes: current fallback route used for final consensus because the scDrugPrio package-example resource universe had zero top-100 overlap with ASGARD across all five states; CLUE `rep_drug_target`/`rep_drug_moa` endpoints returned HTTP 500 on 2026-04-23, so target and MoA annotation now comes from the live `perts` endpoint.

### Additional Output Paths

- `PDOs_outs/Auto_drug_reversal/clue_fallback/Auto_clue_local_rankmatrix_status.csv` — status file for the direct local CMap fallback scoring step.
- `PDOs_outs/Auto_drug_reversal/consensus/Auto_drug_reversal_top5_target_expression.csv` — long-format target-gene expression table used for the final mechanism dot plot across the five finalized PDO states.
####################

####################
### Additional Shell Scripts

- `Auto_run_drug_reversal_local_cmap.sh` — PBS wrapper for `analysis/cell_states/Auto_drug_reversal_local_cmap.R`; uses the dedicated drug-reversal env when present and will propagate `CLUE_API_KEY` or `CLUE_KEY` from the submitted job environment if set.
- `Auto_run_drug_reversal_all_fresh.sh` — orchestrates a full fresh drug-reversal rebuild by forcing `AUTO_DRUG_DEG_MODE=findmarkers`, then submitting dependent ASGARD, scDrugPrio, local CMap fallback, and consensus jobs. Intended for the accurate rerun after DEG or method changes.
####################

####################
### Additional Analysis Scripts

- `analysis/cell_states/Auto_pdo_flot_matched_response.R` — **canonical merged FLOT matched-pair response script**, replacing the separate presentation generation scripts. Produces all original outputs plus three new PDFs:
  - `Auto_pdo_flot_nodeplot_untreated_vs_treated.pdf` — 5-page node-plot PDF: page 1 shows combined (median) untreated vs treated state/hybrid network across all 4 patients; pages 2–5 show per-patient pair node plots.
  - `Auto_pdo_flot_paired_boxplots.pdf` — 3-page paired boxplot PDF: page 1 = state abundance, page 2 = hybrid abundance, page 3 = MP expression; each variable shows untreated (left) vs treated (right) with paired Wilcoxon significance labels.
  - `Auto_pdo_flot_improved_pathway_heatmap.pdf` — improved pathway response heatmap showing mean delta across patients (one column per state), with CC signature score row, lineage MP rows (MP4 Intestinal Metaplasia, MP8 Columnar Progenitor), separate color scales for CC/pathway/lineage, and pairwise significance labels.

### Additional Auto_ Script Dependencies

- `Auto_pdo_flot_matched_response.R`
  Inputs: `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `UCell_scores_filtered.rds`, `Auto_PDO_mp_adj_noreg.rds`, `Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`, `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`; Hallmark gene sets via `msigdbr`
  Env: `dmtcp`
  Notes: replaces separate `Auto_pdo_flot_matched_response.R` + `Auto_pdo_flot_presentation_final.R`. The improved pathway heatmap uses a CC consensus signature (top 50 expressed consensus cell-cycle genes from the Cell_Cycle_Genes.csv reference) and lineage MP expression scored on pseudobulk z-scored logCPM.

### Additional Output Paths

- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_nodeplot_untreated_vs_treated.pdf` — 5-page state/hybrid node-plot network comparing untreated vs FLOT-treated PDOs.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_paired_boxplots.pdf` — 3-page paired boxplot comparing untreated vs treated state abundance, hybrid abundance, and MP expression with Wilcoxon significance.
- `PDOs_outs/Auto_pdo_flot_matched_response/Auto_pdo_flot_improved_pathway_heatmap.pdf` — mean pathway response heatmap with CC signature, lineage MP rows, and significance labels.
####################
