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
Rscript analysis/metaprograms/find_optimal_nmf.R

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

### Analysis Script Governance

All new or substantially updated scripts under `analysis/` must follow the repository map-and-methodology convention:

1. Add a 20-hash `Analysis registry` block at the very start of the script. It must state `Status`, `Script`, `Methodology`, `Map`, exact input files, exact output files or output directories, and whether outputs are downstream inputs or terminal figures only.
2. Reference a methodology file under `analysis/methodology/<matching_analysis_subfolder>/`. The methodology folder structure should mirror `analysis/` where possible.
3. Update `analysis/ANALYSIS_MAP.md` with script status, run order, dependencies, outputs, downstream use, and whether the script is active, terminal, legacy, delete-candidate, or untracked.
4. Update `AGENTS.md` when adding a new script, output path, external reference, recurring technical note, or dependency relationship.
5. Source shared constants/helpers from `analysis/shared/Auto_pdo_analysis_config.R` and `analysis/shared/Auto_pdo_analysis_helpers.R` for new code rather than copying state orders, colors, thresholds, output directory names, metadata column names, cache flags, or logging helpers.
6. Active cleanup renames should use informative names; superseded scripts must use `legacy_` and manual-removal candidates must use `delete_`. Do not delete scripts automatically.
7. If a script is retained only for comparison, mark it as `legacy` and ensure its outputs are not listed as current downstream inputs.

### Output Tiers, Cache/Replot, And Logs

Long-running analysis scripts should write workflow outputs under these subfolders:

- `intermediate/` for heavy RDS/cache/model objects
- `tables/` for final CSV/TSV/XLSX tables
- `figures/` for PDF/PNG figure exports
- `logs/` for lightweight run summaries
- `reports/` for multi-page narrative PDFs or markdown reports

Where possible, scripts should support:

- `PDO_FORCE_REBUILD=1` to ignore cached intermediates and recompute.
- `PDO_REPLOT_ONLY=1` to reuse cached intermediates and regenerate plots/reports only.

Run summaries for heavy scripts should record start/end time, input files, output files, parameters, whether cached objects were reused, and session/package versions when relevant.

### Presentation Plotting

Most downstream plots are used in PowerPoint slides. New figure code must use readable font sizes, legend text, legend keys, row/column labels, point sizes, and line widths at the final exported dimensions. If a visualization is hard to read when inserted into a slide, increase the figure dimensions or split it across pages rather than shrinking labels.

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

### HPC Module Notes
- For Demuxafy/Souporcell reruns on CX3, avoid `module load tools/bioinf` in new PBS wrappers because it can pull stale unavailable dependencies such as `SAMtools/1.15-GCC-11.2.0`. `/usr/bin/singularity` is available after `module purge`; the active Demuxafy wrapper uses that directly with only `module load tools/dev`.
- `BCFtools/1.22-GCC-14.2.0` and `SAMtools/1.22.1-GCC-14.2.0` are under the `tools/prod` module tree (`module spider <module>` reports `tools/prod` as the prerequisite). Load `tools/prod` before these versions.

### ggplot2 + Seurat Compatibility
- **dmtcp env**: ggplot2 was downgraded from 4.0.2 → 3.5.2 because Seurat 5.1.0's `VlnPlot()` crashes internally with ggplot2 4.0.x (`aes_string()` deprecation + `+` operator changes). survminer 0.5.1 and ggpubr remain compatible with ggplot2 3.5.2.
- **gnmf env**: ggplot2 4.0.1 + Seurat 5.4.0 + SeuratObject 5.3.0 — **COMPATIBLE**. Seurat 5.4.0 works with ggplot2 4.0.x (updated internals). No downgrade needed.

### Analysis Directory Structure

```
analysis/
  ANALYSIS_MAP.md  — canonical dependency/run-order/status map for analysis scripts
  shared/          — Auto_pdo_analysis_config.R and Auto_pdo_analysis_helpers.R shared constants, cache controls, output-tier helpers, and logging helpers
  cell_states/     — Auto_legacy_state_hybrid_subtyping_noreg.R (approach B, PDO-adapted), Auto_PDO_state_concordance.R (scRef MPs vs PDO MPs concordance), Auto_pdo_overall_state_proportions.R (overall proportions barplot), Auto_sample_abundance_pdo.R (sample abundance with clinical annotations), Auto_marker_comparison_excel.R (cross-dataset scATLAS vs PDO marker comparison Excel), Auto_five_state_surface_markers.R (PDO-only FACS surface-marker prioritization from five-state marker outputs)
  clinical/        — Auto_clinical_mp_ucell_plots.R, Auto_clinical_variable_plots.R, Auto_survival_clinical_mps.R
  cnv/             — CNV_filter.R, cnv_profile.R, plot_CNV.R (copied from live PDOs)
  enrichment/      — Auto_enrichment_annotation.R, enrichment_extract.R, enrichment_plotting.R, enrich_plot.R, scGSEA.R, wnt_enrich.R
  methodology/     — folder-structured methodology notes mirroring analysis/ where possible
  metaprograms/    — Auto_find_optimal_nmf.R, Auto_extend_nMP_range.R, Auto_update_optimal_mp.R, Auto_mp_correlation_pdo.R, PDO_mp_correlation_crossdata.R, MP_analysis_pdos.R, robust_NMF.R, nmf_plot.R, Find_NMF.R, MP_dist.R, mp_ucell_scoring.R, robust_nmf_scref.R
  plotting/        — heatmap.R
```

Current cleanup map: use `analysis/ANALYSIS_MAP.md` as the authoritative run-order/dependency/status document. The current preferred state route is `Approach B, noreg` via `PDO_states_analysis.R` -> `PDO_unresolved_relabel.R` -> `PDO_finalize_states.R`, with `PDOs_outs/Auto_PDO_final_states.rds` as the preferred downstream state vector. Historical names remain for file safety, but the map records recommended clearer names and legacy candidates.

### Auto_ Script Dependencies

| Script | Input Dependencies | Conda Env |
| :--- | :--- | :--- |
| `Auto_find_optimal_nmf.R` | `Metaprogrammes_Results/geneNMF_metaprograms_nMP_{4..35}.rds` | dmtcp |
| `Auto_extend_nMP_range.R` | `geneNMF_outs.rds` | gnmf |
| `Auto_update_optimal_mp.R` | `MP_outs_default.rds`, `PDOs_merged.rds` | gnmf |
| `Auto_enrichment_annotation.R` | `MP_outs_default.rds` (optimal nMP result) | dmtcp |
| `Auto_mp_correlation_pdo.R` | `PDOs_final.rds`, `MP_outs_default.rds` | dmtcp |
| `Auto_legacy_state_hybrid_subtyping_noreg.R` | `PDOs_merged.rds`, `MP_outs_default.rds`, `UCell_scores_filtered.rds` | dmtcp |
| `Auto_clinical_mp_ucell_plots.R` | `PDOs_merged.rds`, `MP_outs_default.rds`, `UCell_scores_filtered.rds` | dmtcp |
| `Auto_sample_abundance_pdo.R` | `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `UCell_scores_filtered.rds`, Clinical Excel | dmtcp |
| `Auto_clinical_variable_plots.R` | `PDOs_merged.rds` | dmtcp |
| `Auto_survival_clinical_mps.R` | `PDOs_merged.rds`, `UCell_scores_filtered.rds` | dmtcp |
| `PDO_mp_correlation_crossdata.R` | `geneNMF_metaprograms_nMP_13.rds`, `UCell_3CA_MPs.rds`, `scRef/...` | dmtcp |
| `Auto_PDO_state_concordance.R` | `PDOs_final.rds`, `scRef/.../geneNMF_metaprograms_nMP_19.rds`, `Auto_PDO_states_noreg.rds` | dmtcp |
| `Auto_PDO_final_mp_scenic.R` | `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `UCell_scores_filtered.rds`, `UCell_3CA_MPs.rds`, `geneNMF_metaprograms_nMP_13.rds`, `Auto_PDO_unresolved_relabel_mp_coverage.csv`, `New_NMFs.csv` | dmtcp |
| `Auto_marker_comparison_excel.R` | scRef `Auto_six_state_markers_ranked.csv` + `state_specificity.rds`, PDO `Auto_five_state_markers_ranked.csv` + `state_specificity.rds` | dmtcp |
| `Auto_five_state_surface_markers.R` | `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `Auto_five_state_marker_summary.csv`, `Auto_five_state_markers_ranked.csv`, UniProt reviewed-human surface/topology TSV (download/cache), ETH surfaceome Table S3 workbook (download/cache) | dmtcp |
| `Auto_compare_untreated_proportions.R` | `PDOs_all_meta.rds`, `Auto_PDO_final_states.rds` | dmtcp |
| `analysis/cnv/wes_subclone/Auto_run_wes_subclone_sample.sh` | Sarek recalibrated tumour/normal CRAMs, Mutect2 filtered VCFs, `PDOs_outs/Auto_wes_subclone/resources/reference/Homo_sapiens_assembly38.fasta`, `PDOs_outs/Auto_wes_subclone/resources/facets_snps/Auto_ucsc_hg38_snp151Common_biallelic_for_facets.vcf.gz` | local `PDOs_outs/Auto_wes_subclone/conda_env` |
| `analysis/cnv/wes_subclone/Auto_plot_wes_subclone_results.R` | FACETS purity/ploidy and allele-specific segment tables, PyClone-VI input/results, optional conservative Numbat summary | dmtcp |
| `analysis/cnv/wes_subclone/Auto_pyclone_sensitivity.R` | PyClone-VI input tables and local PyClone-VI binary | dmtcp |

### Additional Analysis Scripts

- `analysis/cell_states/Auto_drug_reversal/` — canonical organized drug-reversal workflow folder. It contains the active three-method inhibitor prioritization scripts and wrappers for fresh DEG inputs, ASGARD, scDrugPrio, CLUE/local CMap fallback, reference download/preparation, and visualization generation. Root-level/old copies may remain for file-safety history, but new drug-reversal work should use the folder copies.
- `analysis/cell_states/Auto_pdo_sn_matched_pair_comparison.R` — compares the matched PDO/snRNA-seq pairs `SUR680T3_PDO -> H_post_T1_biopsy` and `SUR791T3_PDO -> L_post_T1_biopsy` using finalized state proportions plus modality-specific top-metaprogram proportions; writes a compact comparison report and CSV summaries to `PDOs_outs/Auto_pdo_sn_matched_pair_comparison/`.
- `analysis/cell_states/Auto_pdo_sn_matched_pair_comparison.R` (update) — now writes a two-page PDF: page 1 focuses on non-cell-cycle state-defining top MPs and renormalized stacked state bars excluding `Unresolved`/`Hybrid`; page 2 keeps the full-view comparison but without subtitles, with centered legends to avoid left-edge clipping.
- `analysis/cell_states/Auto_pdo_flot_matched_response.R` — analyses the four matched untreated/FLOT-treated PDO pairs (`SUR1070`, `SUR1090`, `SUR1072`, `SUR1181`) using finalized state abundance shifts, paired pseudobulk Hallmark-response deltas, heuristic state-fate summaries, and paired edgeR state-specific pseudobulk DE; writes presentation-ready figures plus DEG tables to `PDOs_outs/Auto_pdo_flot_matched_response/`.
- `analysis/cell_states/Auto_pdo_flot_presentation_final.R` — consumes the cached matched-FLOT response outputs and rebuilds the final multi-page presentation PDF, including the pathway matrix, composite response, state abundance, MP change, hybrid-only abundance, recurrent DEG, fate summary, and UMAP support pages.
- `analysis/cell_states/PDO_finalize_states.R` — regenerates `Auto_PDO_final_states.rds` from `unresolved_states/Auto_PDO_unresolved_relabel_states.rds`, merges the stray 3CA respiration/cell-cycle and EMT/protein-maturation labels into the canonical finalized PDO states, and redraws the per-cell heatmap plus state-proportion summaries before the downstream TCGA survival-volcano step.
- `analysis/cell_states/Auto_marker_comparison_excel.R` — canonical cross-dataset marker-comparison workflow for scATLAS vs PDO; includes the combined 3-page marker heatmaps and the top-5 workbook output. The root `Auto_append.R` fragment is redundant.
- `analysis/cell_states/Auto_PDO_scAtlas_scenic_comparison.R` — canonical SCENIC comparison workflow for scATLAS vs PDO; includes RSS gap calculation, the 3-page RSS heatmaps, and the separate top-5 workbook output. The root `Auto_append_scenic.R`, `Auto_fix_rss_gap.R`, and `Auto_update_scenic.R` fragments are redundant.
- `analysis/cell_states/Auto_append_marker_excel.R` — legacy helper fragment for the top-5 workbook layout; its logic is already merged into `Auto_marker_comparison_excel.R`.
- `analysis/cell_states/Auto_compare_untreated_proportions.R` — compares final state proportions (including Unresolved and Hybrid) for SUR1090 and SUR1072 untreated samples; writes a CSV and side-by-side pie charts to `PDOs_outs/Auto_untreated_comparison/`.
- `analysis/cnv/wes_subclone/Auto_*` — reproducible WES tumour-normal subclone workflow for the high-confidence Sarek pairs `PDO_1090_vs_NT_1090` and `PDO_1181_vs_NT_1181`. It prepares a local conda environment, downloads/prepares Broad/GATK GRCh38 reference resources plus UCSC hg38 `snp151Common` for FACETS, validates CRAM `@SQ` SN/LN/M5 compatibility, runs FACETS for purity/ploidy/allele-specific CN, builds/runs PyClone-VI from Mutect2 PASS SNVs, generates visual reliability summaries, and runs a max-cluster-cap PyClone-VI sensitivity check. The lower-confidence pairs `PDO_1070_vs_NT_1070`, `PDO_1072_vs_NT_1072`, `PDO_1121_vs_NT_1121`, and `PDO_1141_vs_NT_1141` are excluded by default.

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
- `PDOs_outs/Auto_wes_subclone/` — WES FACETS -> PyClone-VI workflow outputs, including local software env, downloaded/prepared reference resources, FACETS pileups/segments/purity/ploidy, PyClone-VI inputs/results, and logs for the high-confidence PDO Sarek pairs.
- `PDOs_outs/Auto_wes_subclone/figures/Auto_wes_subclone_visual_summary.pdf` — two-page WES subclone QC visual report with PyClone-VI CCF clusters, raw VAF-vs-CCF, assignment probability, cluster support, FACETS CN segments, and reliability text for `PDO_1090_vs_NT_1090` and `PDO_1181_vs_NT_1181`.
- `PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_subclone_reliability_summary.csv` — per-WES-pair clone-count support metrics and reliability call; `Auto_wes_subclone_cluster_summary.csv` and `Auto_wes_subclone_scRNA_clone_context.csv` provide cluster-level and conservative Numbat context.
- `PDOs_outs/Auto_wes_subclone/tables/visualisation/Auto_wes_subclone_pyclone_sensitivity_summary.csv` — PyClone-VI max-cluster-cap sensitivity table for caps 5, 10, 20, and 40 with three restarts; cluster-level companion table is `Auto_wes_subclone_pyclone_sensitivity_cluster_summary.csv`.

####################

####################
### Additional Analysis Scripts

- `analysis/cell_states/Auto_parse_pdo_mp_scoring_method_trend_check.R` — scores the six Parse trajectory/recovery samples (`T0`, `T1`, `T2`, `T4`, `R4`, `eR4`) with PDO-derived nMP13 metaprograms using three methods: full gene-list UCell, cumulative-weight-filtered UCell, and weighted-rank scoring. Assigns PDO Approach B/noreg states for each scoring method, then evaluates whether lineage MPs/states decrease from `T1` to `T4` and recover at `R4`/`eR4`, while stress MPs/states show the opposite trend. Supports fast-plotting cached execution using saved RDS/CSV outputs when available.

### Additional Auto_ Script Dependencies

- `Auto_parse_pdo_mp_scoring_method_trend_check.R`
  Inputs: Parse per-sample files `/rds/general/project/spatialtranscriptomics/ephemeral/Parse_Pipeline/parse_outs/by_samples/{T0,T1,T2,T4,R4,eR4}/Auto_<sample>_final.rds` with fallback to the `/rds/general/ephemeral/...` Parse path; PDO metaprograms `PDOs_outs/Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`.
  Env: `dmtcp`
  Notes: applies the standard PDO MP filters (`silhouette < 0`, sample coverage < 25%), uses cumulative-weight threshold 0.70 by default, uses Parse count matrices for UCell/rank scoring, and keeps cell-cycle MPs in MP activity plots but excludes them from expected lineage/stress trend scoring. If intermediate scoring files already exist, immediately reuses them to skip expensive computations and regenerate plots.

### Additional Shell Scripts

- `Auto_run_parse_pdo_mp_scoring_method_trend_check.sh` — PBS wrapper for `Auto_parse_pdo_mp_scoring_method_trend_check.R`; uses `dmtcp`, `#PBS -koed`, 4 cores, 128 GB memory, and live log output at `temp/Auto_parse_pdo_mp_scoring_method_trend_check.log`.

### Additional Output Paths

- `PDOs_outs/Auto_parse_pdo_mp_scoring_method_trend_check/Auto_parse_pdo_mp_scores_all_methods.rds` — Parse-cell PDO MP activity matrices for full UCell, cumulative-weight UCell, and weighted-rank scoring.
- `PDOs_outs/Auto_parse_pdo_mp_scoring_method_trend_check/Auto_parse_pdo_state_results_all_methods.rds` and `Auto_parse_pdo_state_assignments_all_methods.csv` — Approach B/noreg state and top-MP assignments for each scoring method.
- `PDOs_outs/Auto_parse_pdo_mp_scoring_method_trend_check/Auto_parse_pdo_mp_activity_sample_summary.csv`, `Auto_parse_pdo_state_abundance_sample_summary.csv`, `Auto_parse_pdo_mp_expected_trend_eval.csv`, `Auto_parse_pdo_state_expected_trend_eval.csv`, and `Auto_parse_pdo_method_expected_trend_summary.csv` — sample-level MP/state summaries and expected-direction scoring tables.
- `PDOs_outs/Auto_parse_pdo_mp_scoring_method_trend_check/Auto_parse_pdo_mp_activity_method_dotplot.pdf/.png`, `Auto_parse_pdo_mp_activity_line_trends.pdf`, `Auto_parse_pdo_state_abundance_method_dotplot.pdf/.png`, `Auto_parse_pdo_state_abundance_line_trends.pdf`, and `Auto_parse_pdo_expected_trend_method_summary.pdf/.png` — comparison figures for choosing the preferred scoring method on Parse data.
####################
### Analysis Cleanup And Governance Additions

- `analysis/ANALYSIS_MAP.md` — canonical map for `analysis/`, including run order, active/legacy/delete-candidate status, input/output dependencies, terminal figure scripts, untracked files that must not be staged, external data requirements, cache/replot policy, and outdated downstream pointers to avoid.
- `analysis/shared/Auto_pdo_analysis_config.R` — central PDO constants for project/output directories, preferred state definition (`Approach B, noreg`), canonical final state vector, state order/colors, MP descriptions, thresholds, metadata column names, plot defaults, external reference paths, and cache environment variable names.
- `analysis/shared/Auto_pdo_analysis_helpers.R` — reusable helper functions for environment-variable parsing, cache policy, output-tier creation, required-file checks, final-state vector name preservation, Seurat assay matrix retrieval, slide-readable ggplot theme defaults, PDF saving, and lightweight run-summary logs.
- `analysis/methodology/README.md` — index and required contents for methodology files.
- `analysis/methodology/shared/shared_config_and_logging_methodology.md` — operational methodology for shared config, output tiers, cache/replot controls, and run summaries.
- `analysis/methodology/cell_states/state_workflows_methodology.md` — operational methodology for the current PDO state-definition route, legacy comparison scripts, and final-state downstream conventions.
- `analysis/methodology/clinical/clinical_association_methodology.md` — operational methodology for clinical association, survival, and final clinical plotting workflows. The single merged script `analysis/clinical/clinical_association_final_figures.R` generates both the final stacked and boxplot clinical figures.
- `analysis/methodology/metaprograms/metaprogram_workflows_methodology.md` — operational methodology for nMP selection, MP filtering, UCell scoring, and MP correlation/cross-dataset comparisons.
- `analysis/methodology/enrichment/enrichment_methodology.md`, `analysis/methodology/cnv/cnv_workflows_methodology.md`, `analysis/methodology/demultiplex/demultiplex_methodology.md`, `analysis/methodology/trajectory/trajectory_methodology.md`, and `analysis/methodology/plotting/plotting_methodology.md` — folder-specific methodology notes for downstream analysis domains.

Current cleanup status:
- Active state-definition files kept with historical names for file safety: `PDO_states_analysis.R`, `PDO_unresolved_relabel.R`, and `PDO_finalize_states.R`. Their headers and `analysis/ANALYSIS_MAP.md` record recommended clearer names.
- Legacy/no-downstream scripts marked in headers and map: `legacy_compare_mp_scoring_state_definition.R`, `legacy_states_scref_pairwise_nodeplot.R`, `legacy_state_hybrid_subtyping_noreg.R`, `legacy_state_hybrid_pairwise_nodeplot_noreg.R`, `legacy_pdo_flot_matched_dge_findmarkers.R`, and `legacy_pdo_flot_matched_survival_and_state_plots.R`.
- `analysis/clinical/clinical_association_final_figures.R` is the canonical, self-contained merged script for generating all final stacked and boxplot clinical figures. Redundant component scripts have been removed.
####################

####################
### Additional Analysis Scripts

- `analysis/cnv/Auto_PDO_numbat_export_inputs.R` — exports per-sample raw-count sparse matrices with raw 10x barcode column names, sample-prefixed cell ID maps, and `Auto_PDO_numbat_manifest.csv` for Numbat. It reuses the PDO velocity/demultiplex CellRanger BAM and QC barcode manifest, excludes `SUR843T3_PDO`, and keeps all outputs under `PDOs_outs/Auto_PDO_numbat/`.
- `analysis/cnv/Auto_PDO_numbat_run_sample.R` — runs Numbat for one sample after `pileup_and_phase.R` has generated allele counts, using the official Numbat container runtime, `ref_hca`, hg38, and `call_clonal_loh=TRUE` for pure malignant PDO samples with no normal-cell compartment.
- `analysis/cnv/Auto_PDO_numbat_concordance_heatmaps.R` — joins Numbat clone calls to existing InferCNA/arm-difference subclone calls, writes per-cell concordance and MP/state composition summaries, and creates one matched-cell PDF page per sample with InferCNA and Numbat heatmaps side by side. `PDO_NUMBAT_CLONE_MODE=conservative` reads the conservative re-cut clone layer and writes separate outputs under `Auto_PDO_numbat/concordance_conservative/`.
- `analysis/cnv/Auto_PDO_numbat_concordance_summary_plots.R` — creates the cohort-level concordance summary PDF from the Numbat-vs-InferCNA tables: metric lollipops, clone-count comparisons, contingency facets, and clone-level state/top-MP composition. `PDO_NUMBAT_CLONE_MODE=conservative` summarizes the conservative concordance tables.
- `analysis/cnv/Auto_PDO_numbat_subclone_mp_heatmap.R` — Numbat-derived subclone MP/state workflow analogous to the InferCNA subclone MP heatmap analysis. It projects Numbat posterior CNV segments onto genome-wide gene bins with uncalled bins set to neutral zero, then plots per-sample Numbat CNV, clone-level MP means/correlations, MP score distributions, final-state abundance, and QC/posterior support. `PDO_NUMBAT_CLONE_MODE=conservative` writes the re-cut clone MP/state outputs under `Auto_PDO_numbat_subclone_mp_conservative/`.
- `analysis/cnv/Auto_PDO_numbat_phylogeny_visualisation.R` — visualizes Numbat's final per-sample phylogenetic trees from `tree_final_<iter>.rds` and `clones_<iter>.rds`, writing one page per sample with the single-cell phylogeny, clone sizes, clone-transition sketch, and clone CNV/event summaries.
- `analysis/cnv/Auto_PDO_numbat_conservative_recut.R` — re-cuts cached Numbat `treeML_<iter>.rds` phylogenies using a conservative `n_cut` strategy, then merges clone branches below `max(20 cells, 3% of cells)` into the best-supported major clone by Numbat posterior probability. This preserves raw Numbat outputs and writes a separate conservative clone layer for 1-4 robust clone groups per sample.

### Additional Auto_ Script Dependencies

- `Auto_PDO_numbat_export_inputs.R`
  Inputs: `PDOs_outs/Auto_velocity_PDO/tables/Auto_pdo_velocity_sample_manifest.csv`, `PDOs_outs/by_samples/<sample>/<sample>.rds`, CellRanger BAMs from `PDOs_outs/Auto_velocity_PDO/cellranger/<sample>/outs/` for Cynthia samples and `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/cellranger/PDOs_{Untreated,Treated}/outs/` for new-batch pools.
  Env: `dmtcp`
- `Auto_PDO_numbat_run_sample.R`
  Inputs: per-sample count RDS/cell map from `Auto_PDO_numbat_export_inputs.R`, per-sample `<sample>_allele_counts.tsv.gz` from Numbat `pileup_and_phase.R`, and the official `pkharchenkolab/numbat-rbase:latest` Singularity image.
  Env: official Numbat container via Singularity/Apptainer.
  Notes: count-matrix cell names are matched to allele-count cells as either raw 10x barcodes or `<sample>_<barcode>` if the pileup output is prefixed.
- `Auto_PDO_numbat_concordance_heatmaps.R`
  Inputs: `PDOs_outs/cnv/Auto_PDO_infercna_target_outs_Carroll_2023.rds`, `PDOs_outs/Auto_PDO_cnv_subclone_mp/Auto_PDO_cnv_subclone_cells.csv`, `PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/Auto_<sample>_numbat_clone_post.csv`, `Auto_<sample>_numbat_joint_post.csv.gz`, optional `Auto_PDO_mp_adj_noreg.rds`.
  Env: `dmtcp`
  Notes: set `PDO_NUMBAT_CLONE_MODE=conservative` to replace only the per-cell Numbat clone labels with `Auto_<sample>_numbat_conservative_clone_post.csv`; Numbat CNV posterior heatmap values still come from the original joint posterior.
- `Auto_PDO_numbat_concordance_summary_plots.R`
  Inputs: `PDOs_outs/Auto_PDO_numbat/concordance/Auto_PDO_numbat_infercna_concordance_summary.csv`, `Auto_PDO_numbat_infercna_contingency.csv`, `Auto_PDO_numbat_clone_state_summary.csv`, and `Auto_PDO_numbat_clone_topmp_summary.csv`.
  Env: `dmtcp`
  Notes: set `PDO_NUMBAT_CLONE_MODE=conservative` to read `PDOs_outs/Auto_PDO_numbat/concordance_conservative/` and write `Auto_PDO_numbat_conservative_concordance_summary_plots.pdf`.
- `Auto_PDO_numbat_subclone_mp_heatmap.R`
  Inputs: `PDOs_outs/Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv`, per-sample normalized Numbat `clone_post` and `joint_post`, `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `UCell_scores_filtered.rds`, `Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`, `Auto_PDO_mp_adj_noreg.rds`, and hg38 gene order.
  Env: `dmtcp`
  Notes: set `PDO_NUMBAT_CLONE_MODE=conservative` to use the robust re-cut clone assignments while reusing original Numbat joint posterior CNV values.
- `Auto_PDO_numbat_phylogeny_visualisation.R`
  Inputs: `PDOs_outs/Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv`, per-sample `tree_final_<iter>.rds`, and per-sample `clones_<iter>.rds` under `PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/`.
  Env: `dmtcp`
- `Auto_PDO_numbat_conservative_recut.R`
  Inputs: `PDOs_outs/Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv`, per-sample `treeML_<iter>.rds`, `geno_<iter>.tsv`, `exp_post_<iter>.tsv`, and `allele_post_<iter>.tsv` under `PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/`.
  Env: official Numbat container via Singularity/Apptainer.
  Notes: default conservative settings are `PDO_NUMBAT_CONSERVATIVE_N_CUT=3`, `PDO_NUMBAT_CONSERVATIVE_MIN_FRAC=0.03`, and `PDO_NUMBAT_CONSERVATIVE_MIN_CELLS=20`; existing tree-cut sweeps are reused unless `PDO_FORCE_REBUILD=1`.

### Additional Shell Scripts

- `analysis/cnv/Auto_prepare_pdo_numbat_container.sh` — PBS wrapper that pulls and validates `pkharchenkolab/numbat-rbase:latest` into `PDOs_outs/Auto_PDO_numbat/Auto_numbat-rbase_latest.sif`; includes `#PBS -koed`.
- `analysis/cnv/Auto_run_pdo_numbat_pileup.sh` — PBS per-sample wrapper for Numbat `pileup_and_phase.R`; uses the sample-specific QC barcode file so multiplexed untreated/treated pool BAMs are processed per biological PDO sample; includes `#PBS -koed`.
- `analysis/cnv/Auto_run_pdo_numbat_sample.sh` — PBS per-sample wrapper for `Auto_PDO_numbat_run_sample.R`; includes `#PBS -koed`.
- `analysis/cnv/Auto_run_pdo_numbat_concordance.sh` — PBS wrapper for the final concordance/heatmap PDF; includes `#PBS -koed`.
- `analysis/cnv/Auto_run_pdo_numbat_summary_plots.sh` — PBS wrapper for the Numbat concordance summary PDF; includes `#PBS -koed`.
- `analysis/cnv/Auto_run_pdo_numbat_subclone_mp.sh` — PBS wrapper for the Numbat subclone MP/state workflow; includes `#PBS -koed`.
- `analysis/cnv/Auto_00_submit_pdo_numbat.sh` — dependency launcher for the full Numbat workflow: prepare container, run per-sample pileup/phasing, run per-sample Numbat, then build concordance heatmaps.

### Additional Output Paths

- `PDOs_outs/Auto_PDO_numbat/Auto_PDO_numbat_manifest.csv` — per-sample BAM, barcode, count, allele-count, and Numbat-output manifest for the 20 included PDO samples.
- `PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/input/Auto_<sample>_counts_raw_barcodes.rds` and `Auto_<sample>_cell_map.csv` — Numbat count matrix and raw-barcode to sample-prefixed cell map.
- `PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/<sample>_allele_counts.tsv.gz` — phased allele counts from Numbat `pileup_and_phase.R`.
- `PDOs_outs/Auto_PDO_numbat/by_samples/<sample>/numbat/Auto_<sample>_numbat_clone_post.csv`, `Auto_<sample>_numbat_joint_post.csv.gz`, and `Auto_<sample>_numbat_segs_consensus.csv` — normalized Numbat outputs for downstream comparison.
- `PDOs_outs/Auto_PDO_numbat/concordance/Auto_PDO_numbat_infercna_matched_heatmaps.pdf` — one page per sample with matched-cell InferCNA and Numbat CNV heatmaps, each split by its own clone/subclone cut and clustered within cut groups.
- `PDOs_outs/Auto_PDO_numbat/concordance/Auto_PDO_numbat_concordance_summary_plots.pdf` — cohort-level visual summary of InferCNA-vs-Numbat concordance, clone counts, contingency structure, and clone-state/top-MP composition.
- `PDOs_outs/Auto_PDO_numbat/concordance/Auto_PDO_numbat_infercna_concordance_summary.csv`, `Auto_PDO_numbat_infercna_contingency.csv`, `Auto_PDO_numbat_clone_state_summary.csv`, `Auto_PDO_numbat_clone_topmp_summary.csv`, and optional `Auto_PDO_numbat_clone_mp_score_summary.csv` — concordance and expression/state comparison tables.
- `PDOs_outs/Auto_PDO_numbat/concordance_conservative/Auto_PDO_numbat_infercna_matched_heatmaps.pdf` — conservative re-cut Numbat clone version of the matched InferCNA-vs-Numbat heatmap report.
- `PDOs_outs/Auto_PDO_numbat/concordance_conservative/Auto_PDO_numbat_conservative_concordance_summary_plots.pdf` — cohort-level conservative concordance summary PDF.
- `PDOs_outs/Auto_PDO_numbat/concordance_conservative/Auto_PDO_numbat_infercna_concordance_summary.csv`, `Auto_PDO_numbat_infercna_contingency.csv`, `Auto_PDO_numbat_clone_state_summary.csv`, `Auto_PDO_numbat_clone_topmp_summary.csv`, and optional `Auto_PDO_numbat_clone_mp_score_summary.csv` — conservative concordance and expression/state comparison tables.
- `PDOs_outs/Auto_PDO_numbat_subclone_mp/Auto_PDO_numbat_subclone_mp_sample_pages.pdf` — one page per sample showing Numbat CNV clone structure with MP and final-state analyses.
- `PDOs_outs/Auto_PDO_numbat_subclone_mp/Auto_PDO_numbat_subclone_cells.csv`, `Auto_PDO_numbat_subclone_summary.csv`, `Auto_PDO_numbat_subclone_mp_tests.csv`, `Auto_PDO_numbat_subclone_mp_subclone_tests.csv`, `Auto_PDO_numbat_subclone_state_tests.csv`, `Auto_PDO_numbat_subclone_compartment_summary.csv`, `Auto_PDO_numbat_subclone_sig_count_summary.csv`, and `Auto_PDO_numbat_subclone_mp_cohort_summary.csv/.pdf` — per-cell clone annotations, sample-level summaries, MP/state association tests, compartment summaries, and cohort-level Numbat subclone MP/state summaries.
- `PDOs_outs/Auto_PDO_numbat_subclone_mp_conservative/Auto_PDO_numbat_subclone_mp_sample_pages.pdf` — conservative re-cut Numbat clone version of the MP/state sample-page report.
- `PDOs_outs/Auto_PDO_numbat_subclone_mp_conservative/Auto_PDO_numbat_subclone_cells.csv`, `Auto_PDO_numbat_subclone_summary.csv`, `Auto_PDO_numbat_subclone_mp_tests.csv`, `Auto_PDO_numbat_subclone_mp_subclone_tests.csv`, `Auto_PDO_numbat_subclone_state_tests.csv`, `Auto_PDO_numbat_subclone_compartment_summary.csv`, `Auto_PDO_numbat_subclone_sig_count_summary.csv`, and `Auto_PDO_numbat_subclone_mp_cohort_summary.csv/.pdf` — conservative per-cell clone annotations, sample-level summaries, MP/state association tests, compartment summaries, and cohort-level clone MP/state summaries.
- `PDOs_outs/Auto_PDO_numbat/phylogeny/Auto_PDO_numbat_phylogenetic_trees.pdf` — one page per sample showing Numbat's final phylogenetic tree visualization and clone-level summaries.
- `PDOs_outs/Auto_PDO_numbat/phylogeny/Auto_PDO_numbat_phylogenetic_tree_summary.csv` — per-sample tree audit table with iteration, vertex/edge/tip counts, clone counts, and largest clone.
- `PDOs_outs/Auto_PDO_numbat/conservative_clones/Auto_PDO_numbat_conservative_phylogenetic_trees.pdf` — one page per sample showing the conservative merged Numbat tree and robust clone-size summary.
- `PDOs_outs/Auto_PDO_numbat/conservative_clones/Auto_PDO_numbat_conservative_clone_summary.csv` — per-sample audit table comparing original Numbat clone counts to conservative re-cut/merged clone counts.
- `PDOs_outs/Auto_PDO_numbat/conservative_clones/Auto_PDO_numbat_tree_cut_sweep.csv` — `n_cut=1..5` sensitivity sweep for each sample.
- `PDOs_outs/Auto_PDO_numbat/conservative_clones/by_samples/<sample>/Auto_<sample>_numbat_conservative_clone_post.csv` — per-cell conservative Numbat clone posterior table with raw re-cut clone columns retained.
####################

####################
### Additional Analysis Scripts

- `analysis/trajectory/Auto_export_pdo_velocity_metadata.R` — exports PDO RNA-velocity metadata with sample-prefixed cell IDs, raw 10x barcodes, Seurat UMAP coordinates, finalized cell states, and pre-unresolved-relabel four-state calls; writes per-sample QC barcode lists and `Auto_pdo_velocity_sample_manifest.csv`.
- `analysis/trajectory/Auto_prepare_pdo_velocity_refs.py` — prepares CellRanger-compatible GRCh38 gene and RepeatMasker GTF files for velocyto under `PDOs_outs/Auto_velocity_PDO/ref/`.
- `analysis/trajectory/Auto_scvelo_pdo_visualise.py` — runs or reloads scVelo independently per PDO sample, uses the existing PDO Seurat UMAP for state-transition calculations, uses a per-sample compressed display embedding for readable plots, writes the core and extended per-sample velocity PDFs, and computes directed four-state velocity alignment edges plus aggregate direction summaries.
- `analysis/trajectory/Auto_pdo_velocity_nodeplots.R` — legacy/deprecated velocity node-plot helper retained for file safety. The active visual summary is now `Auto_pdo_velocity_state_directions.pdf` from `Auto_scvelo_pdo_visualise.py`, so `Auto_run_pdo_scvelo_visualisation.sh` no longer calls this R helper.

### Additional Auto_ Script Dependencies

- `Auto_00_submit_pdo_velocity.sh`
  Inputs: Cynthia FASTQs under `/rds/general/project/spatialtranscriptomics/ephemeral/PDOs/<sample>/`; new-batch demultiplexed CellRanger BAMs under `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/cellranger/PDOs_{Untreated,Treated}/outs/`; PDO `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, and pre-relabel `Auto_PDO_states_noreg.rds`.
  Env: `dmtcp` for metadata export, `velocity` for reference prep/scVelo/velocyto, PBS for CellRanger/filter/velocyto/scVelo.
  Notes: velocity is computed per sample. Cynthia samples require per-sample CellRanger BAM generation from raw FASTQs before filtering and velocyto; new-batch samples reuse the demultiplex rerun pool BAMs and QC barcode lists.

### Additional Shell Scripts

- `analysis/trajectory/Auto_00_submit_pdo_velocity.sh` — full PBS dependency launcher for PDO velocity. It enforces the 46-job throttle, submits per-sample CellRanger for Cynthia samples, filters each sample BAM to QC barcodes, runs per-sample velocyto, then submits dependent scVelo visualisation.
- `analysis/trajectory/Auto_run_pdo_cellranger.sh`, `Auto_filter_sort_pdo_velocity.sh`, `Auto_run_pdo_velocyto.sh`, and `Auto_run_pdo_scvelo_visualisation.sh` — PBS wrappers for the individual PDO velocity stages; each includes `#PBS -koed`.

### Additional Output Paths

- `PDOs_outs/Auto_velocity_PDO/tables/Auto_pdo_velocity_cell_metadata.csv` — per-cell barcode/state/UMAP metadata used to connect velocyto output back to PDO state calls.
- `PDOs_outs/Auto_velocity_PDO/barcodes/<sample>_qc_barcodes.tsv` — per-sample QC barcode lists used to filter BAMs before velocyto.
- `PDOs_outs/Auto_velocity_PDO/cellranger/<sample>/outs/possorted_genome_bam.bam` — CellRanger BAMs generated only for Cynthia-batch single-sample FASTQs.
- `PDOs_outs/Auto_velocity_PDO/coord/<sample>.qc.coord.bam` — per-sample coordinate-sorted BAM filtered to post-QC PDO cells.
- `PDOs_outs/Auto_velocity_PDO/looms/<sample>/` and `PDOs_outs/Auto_velocity_PDO/h5ad/Auto_scvelo_<sample>.h5ad` — per-sample velocyto/scVelo outputs.
- `PDOs_outs/Auto_velocity_PDO/figures/Auto_pdo_velocity_per_sample_visualisations.pdf` — one page per sample with finalized-state UMAP, pre-unresolved-relabel four-state UMAP, and velocity stream only; each page uses one shared finalized-state legend.
- `PDOs_outs/Auto_velocity_PDO/figures/Auto_pdo_velocity_per_sample_visualisations_extended.pdf` — one page per sample with velocity grid, velocity length, and velocity confidence panels.
- `PDOs_outs/Auto_velocity_PDO/figures/Auto_pdo_velocity_state_directions.pdf` — active direction summary: page 1 all PDOs, page 2 Cynthia/new-untreated/new-treated aggregates, then treated/untreated pairs side by side and remaining samples one per page.
- `PDOs_outs/Auto_velocity_PDO/tables/Auto_pdo_velocity_group_state_direction_edges.csv` — aggregate mean/median/fraction-positive four-state direction edges for all PDOs and batch/treatment groups.
- `PDOs_outs/Auto_velocity_PDO/tables/Auto_pdo_velocity_direction_audit_summary.csv` — target-state audit table for checking whether directions are broadly basal-metaplasia-biased or sample/group specific.
- `PDOs_outs/Auto_velocity_PDO/figures/Auto_pdo_velocity_nodeplot_untreated_vs_treated.pdf` — legacy output from the earlier node-plot helper; no longer regenerated by the active scVelo wrapper.
####################

####################
### Additional Analysis Scripts

- `analysis/demultiplex/Auto_01_cellranger_pdo_pool.sh` — PBS CellRanger rerun wrapper for the multiplexed PDO pools (`PDOs_Untreated`, `PDOs_Treated`) using the raw FASTQs under `X204SC25083484-Z01-F001/.../01.RawData`; refuses to overwrite an existing clean output directory.
- `analysis/demultiplex/Auto_02_souporcell_pdo_pool.sh` — PBS Souporcell rerun wrapper that consumes the rerun CellRanger BAM/barcodes, writes a barcode copy instead of gunzipping CellRanger output in place, and runs `k=6` for untreated or `k=4` for treated unless overridden.
- `analysis/demultiplex/Auto_03_reference_and_assign.sh` and `Auto_03_genotyping_save_assign.R` — builds pool-specific normal-WES donor reference genotype VCFs from Strelka `NT_<donor>.strelka.variants.vcf.gz`, exports donor and Souporcell cluster genotypes with BCFtools, then computes `genotyping_save.R`/Demuxafy-style cluster-to-donor Pearson genotype correlations and reciprocal assignment keys. Do not use the older `genotype.sh` overlap-count approach for the current PDO demultiplex rerun.
- `analysis/demultiplex/Auto_04_write_demultiplexed_counts.R` — exports donor-specific count CSVs from a rerun CellRanger matrix plus Souporcell/genotype assignments into the demultiplex rerun output folder, with barcode assignment audit tables.
- `analysis/demultiplex/Auto_05_verify_existing_assignments.R` — diagnostic script for current PDO objects; quantifies existing untreated expression/CNV similarity, UMAP centroid distances, and compares old barcode-to-donor labels against rerun assignments when available.
- `analysis/demultiplex/Auto_compare_SUR1121_SUR1141_wes_cnv.R` — compares donor-level WES CNVkit segment profiles for `PDO_1121_vs_NT_1121` and `PDO_1141_vs_NT_1141` on 1 Mb bins to test whether the suspicious scRNA InferCNA similarity is supported by bulk CNV.
- `analysis/demultiplex/Auto_00_submit_demultiplex_rerun.sh` — login-node submission helper that launches CellRanger, Souporcell, `genotyping_save`-style assignment, count export, and final verification for both multiplexed pools using PBS dependencies.
- `analysis/demultiplex/Auto_README_demultiplex.md` — run-order notes for the organized demultiplex rerun workflow.

### Additional Auto_ Script Dependencies

- `Auto_01_cellranger_pdo_pool.sh`
  Inputs: raw FASTQs in `/rds/general/project/tumourheterogeneity1/live/ITH_sc/X204SC25083484-Z01-F001/X204SC25083484-Z01-F001/01.RawData/<pool>/`; CellRanger `/rds/general/project/tumourheterogeneity1/live/ITH_sc/cellranger-9.0.1/bin/cellranger`; transcriptome `/rds/general/project/tumourheterogeneity1/live/ITH_sc/refdata-gex-GRCh38-2024-A`
  Notes: PBS wrapper requests `select=1:ncpus=16:mem=512gb`, includes `#PBS -koed`, and writes clean rerun outputs under `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/cellranger/`.
- `Auto_02_souporcell_pdo_pool.sh`
  Inputs: rerun CellRanger `possorted_genome_bam.bam` and `filtered_feature_bc_matrix/barcodes.tsv.gz`; genome FASTA `/rds/general/project/tumourheterogeneity1/live/demultiplex/genome.fa`; Demuxafy/Souporcell container from the live demultiplex/multiplexed folders.
  Notes: PBS wrapper requests `select=1:ncpus=18:mem=512gb`, includes `#PBS -koed`, and writes under `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/souporcell/<pool>/`.
- `Auto_03_reference_and_assign.sh`
  Inputs: Souporcell `cluster_genotypes.vcf`; Strelka normal VCFs for `NT_1070`, `NT_1090`, `NT_1072`, `NT_1121`, `NT_1141`, and `NT_1181` as applicable.
  Env: `dmtcp` for the R assignment step; BCFtools module `BCFtools/1.22-GCC-14.2.0`.
  Notes: default donors are `1070,1090,1072,1121,1141,1181` for `PDOs_Untreated` and `1070,1090,1072,1181` for `PDOs_Treated`; per-donor reference VCFs are filtered to biallelic PASS heterozygous SNPs before merge, matching the previous `merged.het.vcf.gz` style used by `genotyping_save.R`.
- `Auto_04_write_demultiplexed_counts.R`
  Inputs: rerun CellRanger matrix, Souporcell `clusters.tsv`, and `Auto_<pool>_cluster_to_donor_key.tsv`.
  Env: `dmtcp`.
- `Auto_05_verify_existing_assignments.R`
  Inputs: current `PDOs_outs/by_samples/*_Untreated_PDO/*.rds`, `PDOs_merged.rds`, optional InferCNA per-sample RDS files, and optional rerun barcode assignment tables.
  Env: `dmtcp`.

### Additional Output Paths

- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/cellranger/<pool>/outs/` — clean CellRanger rerun output for each multiplexed PDO pool.
- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/souporcell/<pool>/` — Souporcell rerun output including `clusters.tsv` and `cluster_genotypes.vcf`.
- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/reference_genotypes/<pool>/Auto_<pool>_donor_reference_snps.vcf.gz` — merged normal-donor genotype reference VCF for assignment.
- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/genotype_assignment/<pool>/Auto_<pool>_cluster_to_donor_key.tsv` and `Auto_<pool>_Genotype_ID_key.txt` — cluster-to-donor and reciprocal donor-to-cluster genotype assignment outputs.
- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_PDO_demultiplex/counts_csv/<pool>/<sample>.csv` — rerun donor-specific count CSVs for QC-pipeline input review before copying into the canonical live count-matrix location.
- `PDOs_outs/Auto_demultiplex_verification/` — current-object and rerun-comparison diagnostics, including the SUR1121/SUR1141 expression/CNV similarity summary.
- `PDOs_outs/Auto_demultiplex_verification/Auto_SUR1121_SUR1141_wes_cnv_summary.csv` and `Auto_SUR1121_SUR1141_wes_cnv_1Mb_profiles.csv` — WES CNVkit 1 Mb profile comparison for SUR1121 vs SUR1141; current summary shows low Pearson correlation despite similar scRNA InferCNA.
- `analysis/demultiplex/Auto_demultiplex_rerun_jobs.tsv` — local submission log with PBS job IDs and dependency edges for the clean demultiplex rerun.
####################

####################
### Additional Analysis Scripts

- `analysis/cnv/Auto_PDO_infercna.R` — PDO-adapted InferCNA workflow based on the Parse `Auto_parse_infercna.R` script. Uses all PDO sample RDS files in `PDOs_outs/by_samples/*_PDO/` except `SUR843T3_PDO`, the Carroll 2023 non-malignant reference, and writes the all-sample CNV-profile heatmap, InferCNA scatter plots, and reusable target/per-sample InferCNA matrices. Supports quick replotting by skipping computation if intermediate files are found.
- `analysis/cnv/Auto_PDO_cnv_subclone_mp_heatmap.R` — PDO-adapted malignant subclone workflow based on scRef `Auto_malignant_subclone_mp_heatmap.R`. Uses the PDO InferCNA target matrix, reconstructs merged PDO cell IDs from `sample + original_cell`, infers per-sample CNA subclones, and summarizes subclone associations with PDO final states and PDO metaprogram scores.
- `analysis/cnv/Auto_PDO_cna_diagnostics_SUR1121_SUR1141.R` — focused audit for the SUR1121/SUR1141 near-duplicate expression-derived CNA result. Compares InferCNA sample means, heatmap-binned InferCNA profiles, raw logCPM pseudobulk expression, and genomic-binned raw expression across untreated samples, and records whether the near-duplicate profile exists before or after InferCNA smoothing.
- `analysis/methodology/cnv/cnv_workflows_methodology.md` — current methodology note for PDO CNV workflows, including InferCNA, CNA subclones, Numbat, output tiers, and cache/replot expectations. The older flat `analysis/methodology/Auto_PDO_cnv_subclone_methodology.md` may remain untracked for file-safety history.

### Additional Auto_ Script Dependencies

- `Auto_PDO_infercna.R`
  Inputs: `PDOs_outs/by_samples/<sample>/<sample>.rds`; Carroll reference `/rds/general/project/tumourheterogeneity1/ephemeral/scRef_Pipeline/ref_outs/Carroll_2023_reference.rds`; gene order `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`
  Env: `dmtcp`
  Notes: PBS wrapper `Auto_pdo_infercna.sh`; excludes `SUR843T3_PDO`; uses Carroll macrophage/endothelial reference groups via `reference$ref`; writes both full reference+target and target-only InferCNA matrices before downstream metrics/plotting. The script now includes quick replot logic and InferCNA scatter visualizations with 99.5th percentile signal capping.
- `Auto_PDO_cnv_subclone_mp_heatmap.R`
  Inputs: `cnv/Auto_PDO_infercna_target_outs_Carroll_2023.rds`, `cnv/Auto_PDO_infercna_target_meta_Carroll_2023.rds`, `PDOs_merged.rds`, `Auto_PDO_final_states.rds`, `Auto_PDO_mp_adj_noreg.rds`, `UCell_scores_filtered.rds`, `Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`; gene order `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`
  Env: `dmtcp`
  Notes: PBS wrapper `Auto_pdo_cnv_subclone_mp.sh`; current schedulable request is `select=1:ncpus=8:mem=64gb`; preserves final-state vector names and maps InferCNA `original_cell` barcodes back to merged PDO cell names before intersecting with state/MP matrices. `percent.mt` is restored from per-sample PDO RDS files because it is absent from `PDOs_merged.rds`; subclone colours are fixed by label; subclone division follows the published Nature-paper style: top 67% CNA-signal genes, Louvain kNN clustering with `k=15`, chromosome-arm CNA calls at ±0.10, dropping provisional clusters below 20 cells or 5% sample fraction, merging same-shape strength-scaled profiles, and otherwise retaining robust arm-level pattern shifts.
- `Auto_PDO_cna_diagnostics_SUR1121_SUR1141.R`
  Inputs: `cnv/Auto_PDO_cnv_sample_mean_profiles_Carroll_2023.rds`, `cnv/Auto_PDO_cnv_heatmap_all_samples_Carroll_2023_input.rds`, `cnv/Auto_PDO_infercna_target_meta_Carroll_2023.rds`, untreated `PDOs_outs/by_samples/<sample>/<sample>.rds`, and gene order `/rds/general/project/spatialtranscriptomics/live/ITH_all/all_samples/hg38_gencode_v27.txt`
  Env: `dmtcp`
  Notes: diagnostic-only script for the SUR1121/SUR1141 question; it does not rerun InferCNA or overwrite existing CNA outputs. It quantifies that the SUR1121/SUR1141 similarity is already present in raw pseudobulk expression and is strongest after genomic binning, so the expression-derived CNA near-duplicate is a method/input-expression artifact rather than barcode duplication.

### Additional Output Paths

- `PDOs_outs/cnv/Auto_PDO_infercna_outs_Carroll_2023.rds` — full InferCNA output matrix including PDO targets and Carroll reference cells.
- `PDOs_outs/cnv/Auto_PDO_infercna_target_outs_Carroll_2023.rds` — PDO target-only InferCNA CNV matrix used for subclone analysis.
- `PDOs_outs/cnv/Auto_PDO_infercna_target_meta_Carroll_2023.rds` and `Auto_PDO_infercna_meta_Carroll_2023.csv` — target/full metadata with sample labels and CNA scatter metrics.
- `PDOs_outs/cnv/Auto_PDO_cnv_heatmap_all_samples_Carroll_2023.pdf` and `_input.rds` — all-PDO-sample binned CNV-profile heatmap and plotting cache.
- `PDOs_outs/cnv/Auto_PDO_cnv_scatter_Carroll_2023.pdf` — InferCNA scatter plots (all-sample and per-sample facets) with signal vs correlation metrics.
- `PDOs_outs/cnv/by_samples/<sample>/Auto_<sample>_infercna_outs_Carroll_2023.rds` — per-sample target InferCNA matrices.
- `PDOs_outs/cnv/Auto_PDO_cna_diagnostics/Auto_SUR1121_SUR1141_cna_diagnostic_summary.csv` — compact metric summary showing SUR1121/SUR1141 InferCNA, binned InferCNA, raw logCPM pseudobulk, and genomic-binned raw-expression correlations plus cell counts and exact-identity checks.
- `PDOs_outs/cnv/Auto_PDO_cna_diagnostics/Auto_infercna_sample_mean_pair_metrics_untreated.csv`, `Auto_infercna_sample_mean_pair_metrics_treated.csv`, `Auto_raw_logcpm_pseudobulk_pair_metrics_untreated.csv`, and `Auto_raw_logcpm_genomic_binned_pair_metrics_untreated.csv` — pairwise profile-similarity diagnostics used to identify SUR1121/SUR1141 as the only near-duplicate untreated expression-derived CNA pair and to confirm treated profiles are distinct.
- `PDOs_outs/Auto_PDO_cnv_subclone_mp/Auto_PDO_cnv_subclone_mp_sample_pages.pdf` — per-sample CNA subclone, MP, state, and QC summary pages.
- `PDOs_outs/Auto_PDO_cnv_subclone_mp/Auto_PDO_cnv_subclone_mp_cohort_summary.pdf` — cohort-level subclone/MP/state/QC summary.
- `PDOs_outs/Auto_PDO_cnv_subclone_mp/Auto_PDO_cnv_subclone_cells.csv`, `Auto_PDO_cnv_subclone_summary.csv`, `Auto_PDO_cnv_subclone_cluster_diagnostics.csv`, `Auto_PDO_cnv_subclone_mp_tests.csv`, `Auto_PDO_cnv_subclone_mp_subclone_tests.csv`, `Auto_PDO_cnv_subclone_state_tests.csv`, `Auto_PDO_cnv_subclone_qc_tests.csv`, `Auto_PDO_cnv_subclone_sig_count_summary.csv`, and `Auto_PDO_cnv_subclone_mp_cohort_summary.csv` — cell-level calls, candidate-cluster diagnostics, and statistical summaries for PDO CNA subclone analysis.
####################

####################
### Additional Analysis Scripts

- `analysis/cell_states/legacy_compare_mp_scoring_state_definition.R` — compares three PDO metaprogram activity/state-definition strategies using the optimal nMP=13 GeneNMF metaprogram object: full gene-list UCell, cumulative-weight-filtered UCell, and weighted-rank activity scoring from `metaprograms.genes.weights`. Reuses the noreg state-definition thresholds from `PDO_states_analysis.R` without unresolved-cell relabelling/finalisation, and writes pairwise activity/state concordance statistics plus PDFs.

### Additional Shell Scripts

- `Auto_run_compare_mp_scoring_state_definition.sh` — PBS wrapper for `legacy_compare_mp_scoring_state_definition.R`; uses `dmtcp`, `#PBS -koed`, 8 cores, 128 GB memory, and live log output at `temp/Auto_compare_mp_scoring_state_definition.log`.
- `Auto_run_compare_mp_scoring_state_definition_4core.sh` — lower-core PBS wrapper for the same analysis when 8-core placement is blocked; uses `dmtcp`, `#PBS -koed`, 4 cores, 96 GB memory, and the same live log output path.

### Additional Auto_ Script Dependencies

- `legacy_compare_mp_scoring_state_definition.R`
  Inputs: `PDOs_merged.rds`, `MP_outs_default.rds` (fallback `Metaprogrammes_Results/geneNMF_metaprograms_nMP_13.rds`), optional `UCell_scores_filtered.rds` for full-UCell audit.
  Env: `dmtcp`
  Notes: excludes `SUR843T3_PDO`; filters MPs with silhouette < 0 and sample coverage < 25%; default cumulative-weight threshold is 0.70 and can be changed with `AUTO_MP_CUM_WEIGHT_THRESHOLD`; rank cap defaults to `AUTO_MP_MAX_RANK=1500`.

### Additional Output Paths

- `PDOs_outs/Auto_compare_mp_scoring_state_definition/Auto_mp_cumulative_weight_gene_counts.csv` and `Auto_mp_gene_weight_selection_long.csv` — MP-level and gene-level cumulative-weight selection diagnostics.
- `PDOs_outs/Auto_compare_mp_scoring_state_definition/Auto_mp_activity_scores_all_methods.rds` — cell-level MP activity matrices for all three scoring methods.
- `PDOs_outs/Auto_compare_mp_scoring_state_definition/Auto_activity_pairwise_scatter.pdf` and `Auto_activity_pairwise_correlation_stats.csv` — three-page pairwise activity scatterplots with per-MP Pearson/Spearman/top-decile concordance statistics.
- `PDOs_outs/Auto_compare_mp_scoring_state_definition/Auto_state_assignment_results_all_methods.rds`, `Auto_PDO_states_<method>.rds`, and `Auto_PDO_top_mp_<method>.rds` — state and top-MP calls for each scoring method.
- `PDOs_outs/Auto_compare_mp_scoring_state_definition/Auto_state_pairwise_concordance_heatmaps.pdf`, `Auto_state_pairwise_concordance_stats.csv`, `Auto_state_pairwise_concordance_counts.csv`, and `Auto_state_composition_by_method.pdf` — pairwise state concordance and composition summaries.
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
- `analysis/methodology/cell_states/Auto_marker_selection_simulation_methodology.md` — persistent methodology note for the marker-selection and qRT-PCR shift simulations.
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
  Notes: runtime can be adjusted with `AUTO_SCDRUGPRIO_CORES`, `AUTO_SCDRUGPRIO_RANDOM_ITERATIONS`, `AUTO_SCDRUGPRIO_PADJ`, `AUTO_SCDRUGPRIO_MIN_ABS_LOGFC`, and `AUTO_SCDRUGPRIO_MAX_DISEASE_GENES`. scDrugPrio ranking must prioritize pharmacological direction before network proximity: `direction_call == "mimicking"` is excluded by default with `AUTO_SCDRUGPRIO_EXCLUDE_MIMICS=1`, and final selected hits require at least one counteracting DEG target with `AUTO_SCDRUGPRIO_REQUIRE_COUNTERACTION=1`. Drugs whose targets have no DEG/direction information should not be treated as significant reversal hits. If network screening succeeds but combined post-processing fails, rerun with `AUTO_SCDRUGPRIO_REUSE_STATE_RESULTS=1` to reuse per-state `Auto_scdrugprio_ranked_<state>.csv` and `Auto_scdrugprio_direction_audit_<state>.csv` files; the script normalizes column types before binding state tables.
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
- `analysis/methodology/cell_states/Auto_drug_reversal_methodology.md` — operational methodology for the ASGARD/scDrugPrio/CLUE consensus reversal workflow.

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

####################
## Nature-Figure Publication Skill (Selective Application)

A `nature-figure` skill is installed at `/rds/general/user/sg3723/home/nature-skills/nature-figure/`. It enforces Nature-journal visual standards. **Agents must exercise judgment to apply this skill primarily to scripts producing final, sharable results.**

### When to Apply (Clever Selection)

Do **not** apply this to every R script. Focus on scripts that synthesize data across samples or produce "Final" visualizations. Prioritize:
- Final cohort-level summary plots (abundance, survival, clinical associations)
- Cross-dataset comparison figures
- Any `Auto_` script producing figures explicitly intended for manuscript inclusion or presentation slides

### How to Apply (R Backend — PDF Priority)

1. **Figure contract**: Define the claim and evidence hierarchy first.
2. **Typography**: Use 6.5pt Arial (Nature standard) via `theme_nature_contract()`.
3. **Export policy (PDF Priority)**: **PDF is the preferred format.** SVG is not required unless requested. Use `grDevices::cairo_pdf()` to ensure font embedding.
   ```r
   save_pub_pdf <- function(plot, filename, width_mm = 183, height_mm = 120) {
     w <- width_mm / 25.4; h <- height_mm / 25.4
     grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
     if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) {
       ComplexHeatmap::draw(plot, merge_legend = TRUE)
     } else {
       print(plot)
     }
     dev.off()
   }
   ```
4. **Color & IA**: Use restrained palettes and follow the **overview → deviation → relationship** information architecture.

### Reference Files
- `~/nature-skills/nature-figure/SKILL.md` — full skill specification
- `~/nature-skills/nature-figure/references/r-workflow.md` — R-specific patterns

### Exceptions
- **Diagnostic/QC scripts**: Step 1-3 pipeline outputs, internal QC heatmaps, and debugging plots should use standard Seurat/ggplot2 defaults to save time.
- **Development/Test scripts**: `delete_*.R` scripts.
####################

####################
### Additional Analysis Scripts

- `analysis/cell_states/Auto_pdo_flot_matched_geneNMF.R` — runs GeneNMF `multiNMF()` on only the eight matched untreated/FLOT-treated PDO samples used in `legacy_pdo_flot_matched_survival_and_state_plots.R` (`SUR1070`, `SUR1072`, `SUR1090`, `SUR1181`; untreated + treated). Writes the matched-sample NMF programme object into one dedicated output folder and does not run nMP selection.
- `analysis/cell_states/Auto_pdo_flot_matched_highres_mp_trend_filter.R` — starts from the matched-sample GeneNMF programme object, sets high-resolution `nMP = round(total NMF programmes / 2)`, runs `getMetaPrograms()`, UCell-scores all cells across the eight matched samples, and retains MPs whose mean and median UCell scores both increase or both decrease in at least three of four treated-vs-untreated pairs. Retained MPs may include single-programme MPs.
- `analysis/cell_states/Auto_pdo_flot_highres_enrichment_annotation.R` — consumes the retained matched-FLOT high-resolution MP genes and trend summary, runs GO BP, Hallmark, 3CA MP, and developmental reference enrichments, and writes a multi-page enrichment PDF plus per-reference/per-column-group PNG heatmaps ordered by paired trend significance.

### Additional Auto_ Script Dependencies

- `Auto_pdo_flot_matched_geneNMF.R`
  Inputs: `PDOs_list_PDOs.rds`
  Env: `gnmf`
  Notes: same `multiNMF()` settings as `geneNMF.R` (`assay="RNA"`, `k=4:9`, `min.exp=0.05`), restricted to the eight matched FLOT PDO samples.
- `Auto_pdo_flot_matched_highres_mp_trend_filter.R`
  Inputs: `Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_matched_geneNMF_outs.rds`; per-sample files `PDOs_outs/by_samples/<sample>/<sample>.rds`; optional 3CA label support from `/rds/general/project/tumourheterogeneity1/live/ITH_sc/PDOs/Count_Matrix/New_NMFs.csv` and cell-cycle genes from `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/Cell_Cycle_Genes.csv`
  Env: `gnmf`
  Notes: trend retention is pairwise and directional: both mean and median treated scores must be greater than untreated for `increase`, or both lower for `decrease`, in at least 3/4 matched pairs. Outputs are sorted by paired mean/median Wilcoxon trend statistics before support count and effect size.
- `Auto_pdo_flot_highres_enrichment_annotation.R`
  Inputs: `Auto_pdo_flot_highres_selected_mp_genes_nMP{k}.rds`, `Auto_pdo_flot_highres_trend_summary_nMP{k}.csv`, 3CA MP reference, MSigDB Hallmark via `msigdbr`, GO via `org.Hs.eg.db`, and developmental references from `/rds/general/project/tumourheterogeneity1/live/EAC_Ref_all/00_merged/developmental/per_stage/*.rds`
  Env: `dmtcp`
  Notes: the `gnmf` environment lacks `clusterProfiler`/`org.Hs.eg.db`; run enrichment after the scoring step in `dmtcp`.

### Additional Shell Scripts

- `Auto_pdo_flot_highres_metaprogram_trends.sh` — PBS wrapper for the matched FLOT high-resolution GeneNMF/trend-filter/enrichment workflow; requests `select=1:ncpus=8:mem=64gb`, uses `#PBS -koed`, runs GeneNMF/UCell in `gnmf`, then switches to `dmtcp` for enrichment annotation.
- `Auto_pdo_flot_highres_enrichment_annotation.sh` — PBS wrapper for rerunning only the matched-FLOT high-resolution enrichment annotation/plotting stage in `dmtcp`; reuses the cached `Auto_pdo_flot_highres_cluster_enrich_nMP{k}.rds` when present.

### Additional Output Paths

- `PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_matched_geneNMF_outs.rds` — matched eight-sample raw GeneNMF programme object.
- `PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_highres_geneNMF_metaprograms_nMP{k}.rds` — high-resolution metaprogram object using `nMP = round(total NMF programmes / 2)`.
- `PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_highres_UCell_scores_nMP{k}.rds` — all-cell UCell score matrix for the high-resolution MPs.
- `PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_highres_trend_summary_nMP{k}.csv` and `Auto_pdo_flot_highres_trend_retained_nMP{k}.csv` — paired untreated-vs-treated trend calls and retained MP subset.
- `PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_highres_activity_boxplots_nMP{k}_selected.pdf`, `Auto_pdo_flot_highres_mean_median_pair_trends_nMP{k}_selected.pdf`, and `Auto_pdo_flot_highres_selected_mean_activity_heatmap_nMP{k}.pdf/.png` — retained-MP visualizations with paired sample spacing and per-patient treated-vs-untreated connections.
- `PDOs_outs/Auto_pdo_flot_highres_metaprogram_trends/Auto_pdo_flot_highres_cluster_enrich_nMP{k}.rds`, `Auto_pdo_flot_highres_enrichment_results_nMP{k}.csv`, `Auto_pdo_flot_highres_enrichment_annotation_nMP{k}.pdf`, and `Auto_pdo_flot_highres_enrich_nMP{k}_<reference>_group*.png` — retained-MP enrichment annotations and heatmaps across all configured reference families.
####################

####################
### Additional Analysis Scripts

- `analysis/trajectory/Auto_PDO_pseudotime_helpers.R` — shared helper library for PDO Monocle3 pseudotime workflows, adapted from the Parse/scRef trajectory scripts. It uses `PDOs_merged.rds` plus the pre-unresolved-relabel four-state vector `Auto_PDO_states_noreg.rds`, processes each sample independently without batch correction, and defaults the trajectory root to `Basal to Intest. Meta` (`AUTO_PDO_ROOT_STATE` can override this).
- `analysis/trajectory/Auto_PDO_pseudotime_samples.R` — builds and caches per-sample Monocle3 trajectories for valid PDO samples, writes per-cell pseudotime metadata, and exports a combined state/pseudotime PDF.
- `analysis/trajectory/Auto_PDO_pseudotime_linear_plot.R` — reuses cached per-sample trajectories to generate Parse/scRef-style four-panel trajectory reports with principal-graph projections, UMAP pseudotime, root labels, and count-weighted state ridges.
- `analysis/trajectory/Auto_PDO_pseudotime_state_distance_matrix.R` — computes per-sample four-state distances, then summarizes directed pseudotime, principal-graph geodesic, and UMAP centroid distances across samples.

### Additional Auto_ Script Dependencies

- `Auto_PDO_pseudotime_samples.R`
  Inputs: `PDOs_merged.rds`, `Auto_PDO_states_noreg.rds`
  Env: `dmtcp`
  Notes: uses the four pre-relabel states only (`Classic Proliferative`, `Basal to Intest. Meta`, `SMG-like Metaplasia`, `Stress-adaptive`); default inclusion thresholds are 80 total primary-state cells, at least 2 states with 20 cells, and at least 20 root-state cells.
- `Auto_PDO_pseudotime_linear_plot.R`
  Inputs: cached assets from `PDOs_outs/Auto_PDO_pseudotime_pre_relabel/sample_trajectory_assets/`; rebuilds missing assets through the helper if needed.
  Env: `dmtcp`
- `Auto_PDO_pseudotime_state_distance_matrix.R`
  Inputs: cached per-sample Monocle3 `cds`, pseudotime, projection, and metadata assets from `Auto_PDO_pseudotime_samples.R`
  Env: `dmtcp`

### Additional Shell Scripts

- `analysis/trajectory/Auto_run_PDO_pseudotime.sh` — PBS wrapper for the full PDO pseudotime workflow; requests `select=1:ncpus=8:mem=96gb`, includes `#PBS -koed`, activates `dmtcp`, and runs the sample, linear-report, and state-distance scripts in order.

### Additional Output Paths

- `PDOs_outs/Auto_PDO_pseudotime_pre_relabel/Auto_PDO_pseudotime_combined.pdf` — per-sample Monocle3 state and pseudotime plots.
- `PDOs_outs/Auto_PDO_pseudotime_pre_relabel/Auto_PDO_pseudotime_linear_reports.pdf` — multi-page Parse/scRef-style trajectory report, one page per valid PDO sample.
- `PDOs_outs/Auto_PDO_pseudotime_pre_relabel/Auto_PDO_pseudotime_metadata.csv` and `Auto_PDO_pseudotime_summary.csv` — per-cell and per-sample/per-state pseudotime summaries.
- `PDOs_outs/Auto_PDO_pseudotime_pre_relabel/sample_trajectory_assets/Auto_PDO_<sample>_cds.rds`, `_pseudotime.rds`, `_metadata.csv`, `_projections.csv` — cached per-sample trajectory assets.
- `PDOs_outs/Auto_PDO_pseudotime_pre_relabel/state_distance_pseudotime/Auto_PDO_state_distance_summary.csv` and `Auto_PDO_state_distance_matrices.rds` — aggregated four-state distance summaries and matrices.
- `PDOs_outs/Auto_PDO_pseudotime_pre_relabel/state_distance_pseudotime/Auto_PDO_state_distance_method_comparison_heatmap.pdf` and `Auto_PDO_state_distance_nodeplot.pdf` — four-state distance visualizations.
####################
