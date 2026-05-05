# PDO Drug-Reversal Methodology

This document describes how the current PDO drug-reversal workflow is implemented in code.

It is written against the current files:

- `analysis/cell_states/Auto_drug_reversal_inputs.R`
- `analysis/cell_states/Auto_drug_reversal_asgard.R`
- `analysis/cell_states/Auto_drug_reversal_scdrugprio.R`
- `analysis/cell_states/Auto_drug_reversal_local_cmap.R`
- `analysis/cell_states/Auto_drug_reversal_consensus_visuals.R`
- `analysis/cell_states/Auto_drug_reversal_method_visuals.R`
- `analysis/cell_states/Auto_drug_reversal_predicted_reversion_visuals.R`
- `analysis/cell_states/Auto_prepare_asgard_reference.R`
- `Auto_download_asgard_l1000_reference.sh`
- `Auto_build_asgard_reference.sh`

The description below is intentionally operational. It documents what the pipeline is doing, which files it reads and writes, how each score is created, and how each output should be interpreted.

## 1. Biological Aim

The goal of this workflow is to identify compounds that reverse the transcriptional program of each finalized malignant PDO state.

The strategy is not to predict generic cytotoxicity. Instead, it asks a more specific question:

> For a given PDO malignant state, which perturbagens push the state-defining expression program in the opposite direction?

The workflow therefore uses transcriptomic reversal rather than cell-line viability prediction as the primary evidence type.

The implemented design has three layers:

1. Build fresh state-vs-rest differential expression signatures from the finalized five-state PDO object.
2. Score those signatures against large perturbational reference resources.
3. Prioritize compounds that are consistently ranked across independent reversal frameworks, then inspect whether their known target genes are actually expressed in the tumour cells.

## 2. Fixed State Definitions

The workflow is hard-coded to the five finalized PDO malignant states:

- `Classic Proliferative`
- `Basal to Intest. Meta`
- `SMG-like Metaplasia`
- `Stress-adaptive`
- `3CA_EMT_and_Protein_maturation`

These are stored in `state_order` across the drug-reversal scripts and define:

- the DE contrasts
- the exported signature names
- the plotting order
- the per-state consensus outputs

The workflow excludes:

- `SUR843T3_PDO`
- `Hybrid`
- `Unresolved`
- any additional non-finalized labels outside the five-state set

This is deliberate. The inhibitor search is restricted to the finalized malignant states only.

For display, the current canonical order is `Classic Proliferative`, `Basal to Intest. Meta`, `SMG-like Metaplasia`, `Stress-adaptive`, and `3CA_EMT_and_Protein_maturation`. If a scATLAS-style state set includes `Immune Infiltrating`, it is inserted after `Stress-adaptive` and before the 3CA EMT/protein-maturation state.

## 3. Data Sources

### 3.1 Core PDO inputs

The workflow reads:

- `PDOs_outs/PDOs_merged.rds`
- `PDOs_outs/Auto_PDO_final_states.rds`

It will preferentially reuse a cached five-state object if present:

- `PDOs_outs/Auto_drug_reversal/cache/Auto_drug_reversal_state5.rds`

If that cache is absent, or rebuilding is forced, it reconstructs the five-state object from the merged PDO Seurat object plus the finalized state vector.

### 3.2 Optional cached five-state object from the marker workflow

If available, the input script can also reuse:

- `PDOs_outs/Auto_five_state_markers/cache/pdos_state5_embedded.rds`

This speeds up input generation because the five-state subset and its normalized assay are already present.

### 3.3 External perturbation resources

Three external perturbational resources are used:

1. ASGARD tissue-specific LINCS rank-matrix references.
2. scDrugPrio network resources:
   - protein-protein interaction network
   - drug-target mapping table
   - optional pharmacological action table
3. CLUE perturbagen metadata through the live `perts` endpoint for compound target and MoA annotation.

## 4. Step 1: Building the Five-State Drug-Reversal Object

This is implemented in `Auto_drug_reversal_inputs.R`.

### 4.1 Object construction

If rebuilding is needed, the script:

1. loads `PDOs_merged.rds`
2. loads `Auto_PDO_final_states.rds`
3. intersects cells present in both
4. keeps only cells whose state label is one of the five finalized states
5. removes `SUR843T3_PDO`

The retained cells are turned into a dedicated Seurat object for the drug-reversal workflow.

### 4.2 Metadata retained

The script attaches:

- `orig.ident`
- `Batch` if present
- derived `batch`
- `state`

The derived `batch` field is set by sample naming pattern:

- samples containing `_Treated_` or `_Untreated_` -> `New_batch`
- all others -> `Cynthia_batch`

This batch label is not currently used as a DE covariate in the implemented code. It is exported so downstream wrappers and QC summaries can track sample origin.

### 4.3 Assay used

The workflow uses the `RNA` assay.

When rebuilding from the merged object:

1. the `counts` layer is extracted
2. genes detected in fewer than 10 cells are removed
3. a new Seurat object is created from the filtered count matrix
4. `NormalizeData()` is run
5. the cell identities are set to the five-state factor

This means the DEG step is run on a state-restricted RNA assay rather than on the full merged multi-label object.

### 4.4 Summary exports from object preparation

Before DE is run, the script writes:

- `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_state_sample_cell_counts.csv`
- `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_metadata.csv`

These files record:

- how many cells each sample contributes to each state
- which state each retained cell belongs to

Interpretation:

- use `Auto_drug_reversal_state_sample_cell_counts.csv` to check whether a state is dominated by a small number of samples
- use `Auto_drug_reversal_metadata.csv` to confirm exactly which cells entered the reversal workflow

## 5. Step 2: Fresh State-vs-Rest Differential Expression

This is the key input-generation step. It creates the disease signatures that every downstream drug-reversal method uses.

### 5.1 Contrast design

For each of the five states, the script compares:

- `ident.1 = target state`
- `ident.2 = all other finalized malignant states`

This is a state-vs-rest comparison inside the malignant compartment only.

It is not:

- state vs all cells in the full dataset
- state vs matched normal
- pairwise state vs one other state at a time

### 5.2 DE method

The implemented default fresh-DE path is:

- `FindMarkers()`
- `test.use = "wilcox"`
- `logfc.threshold = 0`
- `min.pct = 0.01` unless overridden by `AUTO_DRUG_DEG_MIN_PCT`
- `only.pos = FALSE`

Because `test.use = "wilcox"` is used, Seurat will leverage `presto` when the package is installed in the active environment. That is the intended fast path for the fresh DEG rerun.

### 5.3 Why `logfc.threshold = 0`

The drug-reversal step needs a complete ranked signature, not only highly filtered positive markers.

For that reason the script does not pre-truncate the DE table at the testing stage. It keeps:

- strongly positive genes
- strongly negative genes
- weakly shifted genes

The ranking and top-signature extraction happen later.

### 5.4 DEG cache and force-rerun controls

The full DEG table is written to:

- `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_degs_all_states.csv.gz`

If that file already exists, the script will reuse it unless:

- `AUTO_FORCE_DRUG_DEGS=1`

This is how the accurate rerun is triggered after method upgrades.

### 5.5 Fallback DEG mode

If `AUTO_DRUG_DEG_MODE=global`, the script does not run `FindMarkers()`.

Instead it uses:

- `PDOs_outs/Auto_five_state_markers/Auto_five_state_global_marker_screen.csv.gz`

This fallback was added so the reversal workflow could still be developed when the fresh pooled Wilcoxon step was too slow in the initial environment.

Important interpretation point:

- `AUTO_DRUG_DEG_MODE=findmarkers` = statistical fresh state-vs-rest DE
- `AUTO_DRUG_DEG_MODE=global` = descriptive fallback based on an existing global marker screen

The accurate production run should use the fresh `findmarkers` path.

## 6. Step 3: Deriving Ranked Signatures From the Full DEG Table

Once the full DEG table exists, the script derives compact signature files for each downstream framework.

### 6.1 Internal ranking columns

The script computes:

- `p_rank_value = p_val_adj`, using `1` if `p_val_adj` is missing
- `abs_logfc = abs(avg_log2FC)`

These are used for ordering the DEG table.

### 6.2 Upregulated signature selection

For each state, upregulated genes are:

1. restricted to `avg_log2FC > 0`
2. ordered by:
   - ascending adjusted p-value
   - descending log fold-change
   - descending `pct.1`
3. truncated to the top `AUTO_DRUG_SIGNATURE_TOP_N`, default `150`

### 6.3 Downregulated signature selection

For each state, downregulated genes are:

1. restricted to `avg_log2FC < 0`
2. ordered by:
   - ascending adjusted p-value
   - ascending log fold-change (most negative first)
   - descending `pct.2`
3. truncated to the top `AUTO_DRUG_SIGNATURE_TOP_N`, default `150`

### 6.4 Why top 150 up and top 150 down

This size is a compromise between:

- retaining enough biology to represent the state
- avoiding overly broad signatures that become noisy in perturbational matching

The same top-150 format is used consistently for:

- CLUE-style queries
- local CMap fallback scoring
- concise state-level presentation summaries

### 6.5 Signature export

The combined per-state top signatures are written to:

- `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_signature_top150.csv`

This file is the central compact signature file for the workflow.

Interpretation:

- each row is one gene assigned to one state and one direction
- `direction = up` means the gene is elevated in the target state
- `direction = down` means the gene is depleted in the target state
- `signature_rank` is the within-direction ordering used for export

## 7. Step 4: Wrapper-Specific Input Exports

The fresh DEG table is reformatted into multiple downstream input styles.

### 7.1 ASGARD inputs

The script builds a named gene-list object for ASGARD:

- `PDOs_outs/Auto_drug_reversal/asgard_inputs/Auto_asgard_gene_list.rds`

For each state the ASGARD table contains:

- row names = genes
- `score = avg_log2FC`
- `adj.P.Val`
- `P.Value`

State-specific tabular versions are also written:

- `PDOs_outs/Auto_drug_reversal/asgard_inputs/Auto_asgard_deg_<state>.txt`

### 7.2 scDrugPrio inputs

The script also writes per-state DEG tables for scDrugPrio:

- `PDOs_outs/Auto_drug_reversal/scdrugprio_inputs/Auto_scdrugprio_deg_<state>.txt`

These include:

- `gene`
- `p_val`
- `avg_logFC`
- `pct.1`
- `pct.2`
- `p_val_adj`

### 7.3 CLUE/CMap GMT inputs

Four GMT files are created:

- `Auto_clue_up_symbols.gmt`
- `Auto_clue_down_symbols.gmt`
- `Auto_clue_up_entrez.gmt`
- `Auto_clue_down_entrez.gmt`

The symbol-to-Entrez mapping is also exported:

- `Auto_clue_symbol_to_entrez.csv`

Interpretation:

- symbol GMT files are human-readable and useful for QA
- Entrez GMT files are the submission-ready form for CLUE/L1000-style APIs that expect Entrez identifiers

### 7.4 Sparse matrix export

If matrix export is enabled, the script writes:

- `matrix/Auto_drug_reversal_counts.mtx`
- `matrix/Auto_drug_reversal_features.tsv`
- `matrix/Auto_drug_reversal_barcodes.tsv`

These are intended for wrappers or external tools that need a sparse matrix plus feature/barcode lists.

## 8. Step 5: ASGARD Reference Preparation

ASGARD cannot run from the PDO DEG list alone. It also needs tissue-specific LINCS reference files.

### 8.1 Raw download

`Auto_download_asgard_l1000_reference.sh` downloads the required GEO LINCS resources into:

- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000/`

The raw staging includes:

- `GSE70138`
- `GSE92742`
- associated cell, gene, and signature metadata

### 8.2 Reference build

`Auto_prepare_asgard_reference.R` runs `Asgard::PrepareReference()` on the uncompressed staging files and writes tissue-specific files under:

- `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000/DrugReference/`

The resulting paths are recorded in:

- `PDOs_outs/Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.csv`
- `PDOs_outs/Auto_drug_reversal/asgard_reference/Auto_asgard_reference_paths.sh`

### 8.3 Tissue choice

The implemented default tissue is:

- `stomach`

This is a proxy choice, not a claim that stomach is biologically identical to oesophageal adenocarcinoma.

It is used because the ASGARD/LINCS reference metadata available in this workflow did not expose an oesophagus primary site.

Interpretation:

- ASGARD is being run on the best available upper-GI epithelial proxy in the available reference resource
- this should be stated when presenting results

## 9. Step 6: ASGARD Drug-Reversal Scoring

This is implemented in `Auto_drug_reversal_asgard.R`.

### 9.1 Inputs

The script loads:

- `Auto_asgard_gene_list.rds`
- either a prebuilt ASGARD reference RDS, or
- `AUTO_ASGARD_DRUG_RESPONSE`
- `AUTO_ASGARD_GENE_INFO`
- `AUTO_ASGARD_DRUG_INFO`

### 9.2 Core call

ASGARD is run through:

- `Asgard::GetDrugRef()` if needed to build the in-memory reference
- `Asgard::GetDrug()`

with:

- `repurposing.unit = "drug"`
- `connectivity = "negative"`
- `drug.type = "all"` by default

### 9.3 Meaning of negative connectivity

Negative connectivity means ASGARD is looking for perturbagens whose transcriptional effects oppose the state-defining DEG pattern.

Operationally, this means:

- genes up in the target state should be driven down by the drug
- genes down in the target state should be driven up by the drug

### 9.4 ASGARD output standardization

Because ASGARD can return differently structured result tables depending on version and internals, the wrapper standardizes the outputs into:

- `state`
- `pipeline`
- `drug`
- `rank`
- `score`
- `p_value`
- `fdr`
- `moa`
- `target_genes`
- `rank_direction`

These are written to:

- `PDOs_outs/Auto_drug_reversal/asgard/Auto_asgard_ranked_drugs.csv`

### 9.5 How to interpret `Auto_asgard_ranked_drugs.csv`

- one row = one drug-state prediction
- lower `rank` = stronger ASGARD priority
- smaller `fdr` or `p_value` = stronger statistical support if present
- `target_genes` and `moa` are optional annotations carried through when ASGARD provides them

## 10. Step 7: scDrugPrio Network Reversal

This is implemented in `Auto_drug_reversal_scdrugprio.R`.

### 10.1 What this module is doing

The workflow uses the scDrugPrio transcriptomic/network prioritization path only.

It does not run:

- the scDrug cell-line cytotoxicity module
- dose-response viability prediction

### 10.2 Resource loading

The script expects:

- `AUTO_SCDRUGPRIO_PPI`
- `AUTO_SCDRUGPRIO_DRUG_TARGETS`
- optional `AUTO_SCDRUGPRIO_PHARMA_EFFECT`

If explicit resources are absent, it falls back to the package example resources:

- `lit_ppi`
- `drug_bank_example_data`

This fallback is useful for testing the code path, but it is not equivalent to a full curated production drug-target universe.

### 10.3 Function-source fallback

If the `scDrugPrio` package is not installed, the wrapper downloads and sources the required raw R functions from the upstream repository into:

- `PDOs_outs/Auto_drug_reversal/resources/scdrugprio_source/`

This avoids blocking the workflow on package installation problems.

### 10.4 Identifier harmonization

The wrapper checks whether the PPI network uses:

- gene symbols, or
- Entrez IDs

If the network is Entrez-based, the DEG genes and drug targets are mapped from symbols to Entrez using:

- `AnnotationDbi`
- `org.Hs.eg.db`

### 10.5 Disease gene selection for each state

For each state, the wrapper starts from the per-state DEG table and selects disease genes as:

- `p_val_adj <= 0.05`
- `abs(avg_logFC) >= 0.10`

If that leaves fewer than 10 genes, it falls back to the top 300 genes ranked by:

- adjusted p-value
- absolute log fold-change

### 10.6 Network screen

The core scDrugPrio function computes drug-target to disease-gene proximity on the PPI network using:

- observed closest-distance metrics
- random degree/bin-adjusted comparisons
- repeated randomization, default `AUTO_SCDRUGPRIO_RANDOM_ITERATIONS = 1000`

The wrapper stores the upstream output files in a per-state directory under:

- `PDOs_outs/Auto_drug_reversal/scdrugprio/<state>/`

### 10.7 Drug-level ranking logic

After the network screen, the wrapper keeps rows with:

- finite `p_value`
- `p_value < 0.05`
- finite `dc`
- `dc < 1`

If none pass, it falls back to the top 500 finite rows ranked by:

- `p_value`
- `dc`

### 10.8 Pharmacological directionality

If a pharmacological action table is available, the wrapper labels targets as:

- counteracting
- mimicking

based on:

- whether the target state DEG is up or down
- whether the drug action implies activation or inhibition

It then computes:

- `counteracting_targets`
- `mimicking_targets`
- `n_counteracting_targets`
- `n_mimicking_targets`
- `counteract_fraction`

### 10.9 Final scDrugPrio score

The wrapper defines:

- `score = -zc + counteract_fraction`

This means higher priority is given to drugs that have:

- stronger favorable network proximity (`zc` more negative)
- greater fraction of targets acting in the opposite direction of the disease signature

### 10.10 scDrugPrio outputs

Per-state ranked files:

- `PDOs_outs/Auto_drug_reversal/scdrugprio/<state>/Auto_scdrugprio_ranked_<state>.csv`

Combined file:

- `PDOs_outs/Auto_drug_reversal/scdrugprio/Auto_scdrugprio_ranked_drugs.csv`

Status file:

- `PDOs_outs/Auto_drug_reversal/scdrugprio/Auto_scdrugprio_status.csv`

### 10.11 How to interpret `Auto_scdrugprio_ranked_drugs.csv`

- lower `rank` = higher scDrugPrio priority
- more favorable `p_value` and lower `dc` indicate stronger network closeness
- more negative `zc` indicates stronger enrichment relative to random expectation
- larger `counteract_fraction` means more of the mapped targets oppose the disease-program direction

Important caveat:

- if the run uses the package example resources, treat the output as a code-path validation or rough prioritization, not as the final consensus evidence layer

## 11. Step 8: Direct Local CMap Fallback

This is implemented in `Auto_drug_reversal_local_cmap.R`.

This is the wrapper-bypass fallback route. It directly scores the PDO signatures against the ASGARD tissue-specific LINCS rank matrix.

### 11.1 Inputs

The script loads:

- `Auto_drug_reversal_signature_top150.csv`
- `Auto_asgard_reference_paths.csv`

and from the ASGARD reference paths:

- `stomach_rankMatrix.txt`
- `stomach_gene_info.txt`
- `stomach_drug_info.txt`

### 11.2 Gene-level rank matrix reconstruction

The rank matrix is stored by probe.

The script:

1. joins the rank matrix to probe-to-gene annotation
2. removes missing gene symbols
3. averages multiple probes mapping to the same gene

This yields a gene-by-instance rank matrix.

### 11.3 Rank normalization

The rank matrix is converted to a normalized scale:

- `norm = (rank - 1) / (n_genes - 1)`

This rescales every perturbation instance so ranks are comparable on a 0 to 1 scale.

Interpretation:

- low normalized rank = gene is driven toward the bottom of the ranked profile
- high normalized rank = gene is driven toward the top of the ranked profile

### 11.4 Reversal score formula

For each state and each perturbation instance, the script computes:

- mean normalized rank of the state up-signature
- mean normalized rank of the state down-signature

and defines:

- `instance_score = mean(norm_rank(up_genes)) - mean(norm_rank(down_genes))`

Interpretation:

- larger positive score means up-signature genes are pushed high while down-signature genes are pushed low
- because this is used as a reversal score in the current implementation, top-ranked compounds are those with the strongest separation under this ranking convention

At the compound level, the script averages across all matching perturbation instances:

- `score = mean(instance_score)`

It also stores:

- `best_instance_score`
- `n_instances`

### 11.5 Compound naming

The script maps `instance_id` to compound names using `stomach_drug_info.txt`.

When `cmap_name` and `catalog_name` are both available, the catalog suffix is stripped so the exported drug label is closer to the conventional perturbagen name.

### 11.6 CLUE metadata annotation

If `CLUE_API_KEY` or `CLUE_KEY` is available, the script annotates top compounds through:

- `https://api.clue.io/api/perts`

It requests:

- `pert_iname`
- `target`
- `moa`

These are used to populate:

- `target_genes`
- `moa`

This `perts` route is the current live endpoint used by the workflow for metadata annotation.

### 11.7 Local CMap outputs

The ranked output is:

- `PDOs_outs/Auto_drug_reversal/clue_fallback/Auto_clue_ranked_drugs.csv`

Status is recorded in:

- `PDOs_outs/Auto_drug_reversal/clue_fallback/Auto_clue_local_rankmatrix_status.csv`

### 11.8 How to interpret `Auto_clue_ranked_drugs.csv`

- one row = one drug-state prediction from the direct local LINCS fallback
- lower `rank` = higher local fallback priority
- higher `score` = stronger signature separation under the implemented local rank-matrix scoring rule
- `n_instances` tells you how many perturbation instances were aggregated into the drug-level score
- `target_genes` and `moa` come from CLUE perturbagen metadata when available

## 12. Step 9: Consensus Prioritization

This is implemented in `Auto_drug_reversal_consensus_visuals.R`.

### 12.1 Ranking sources considered

The consensus script reads:

- `Auto_asgard_ranked_drugs.csv`
- `Auto_scdrugprio_ranked_drugs.csv`
- optional `Auto_clue_ranked_drugs.csv`

### 12.2 Standardization

Each input ranking is standardized into:

- `state`
- `drug`
- `rank`
- `score`
- `p_value`
- `fdr`
- `target_genes`
- `moa`
- `drug_key = tolower(trimws(drug))`

This harmonization is required because the three pipelines do not emit identical column names or structures.

### 12.3 Default secondary pipeline choice

The script prefers:

- ASGARD as the primary ranking
- scDrugPrio as the secondary ranking

### 12.4 Automatic fallback replacement

If scDrugPrio exists but its top-100 overlap with ASGARD is zero across states, and the local CLUE/CMap fallback exists, the script automatically switches the secondary ranking to:

- `CLUE_FALLBACK_LOCAL`

This is how the workflow avoids producing empty consensus figures when the available scDrugPrio resource universe is not commensurate with the ASGARD/LINCS space.

### 12.5 Top-100 overlap summary

For each state the script computes:

- `asgard_top100_n`
- `secondary_top100_n`
- `overlap_n`
- `fallback_used`
- `secondary_pipeline`

This is written to:

- `PDOs_outs/Auto_drug_reversal/consensus/Auto_drug_reversal_top100_overlap_summary.csv`

Interpretation:

- `overlap_n` is the simplest high-confidence consensus count
- `fallback_used = TRUE` means the local CMap fallback replaced scDrugPrio in the final figures

### 12.6 Consensus table

The final consensus table is the inner join of:

- ASGARD top 100 per state
- secondary pipeline top 100 per state

For each matched compound, the script computes:

- `mean_rank = (asgard_rank + scdrugprio_rank) / 2`
- `rank_sum = asgard_rank + scdrugprio_rank`

The rows are then ordered by:

1. `rank_sum`
2. `mean_rank`
3. `drug`

and assigned:

- `consensus_rank`

This table is written to:

- `PDOs_outs/Auto_drug_reversal/consensus/Auto_drug_reversal_consensus_drugs.csv`

### 12.7 State-specificity margin

For each consensus compound, the script also computes how well the prioritization is confined to the target state rather than appearing equally high in off-target states.

It does this by finding the best off-state rank in:

- ASGARD
- the secondary pipeline

and then computing:

- `best_off_state_rank`
- `state_specificity_margin = best_off_state_rank - mean_rank`

Interpretation:

- positive margin = the compound ranks better in the target state than in any off-target state
- negative margin = the compound is ranked at least as strongly for some other state and is therefore less state-specific

## 13. Step 10: Visualizations

The workflow creates three presentation-ready visual outputs for each five-state analysis.

## 13.1 Venn Diagram: top-100 overlap

Output:

- `Auto_drug_reversal_venn_top100.pdf`

Generation:

For each state, the script takes:

- the top 100 ASGARD compounds
- the top 100 secondary-pipeline compounds

and plots a two-set Venn diagram showing:

- size of the ASGARD set
- size of the secondary set
- size of the intersection

Interpretation:

- the overlap count is the immediate high-confidence candidate pool
- a larger intersection suggests more agreement between reversal frameworks
- a very small intersection suggests either a highly state-specific signal or a mismatch between perturbational universes/resources

## 13.2 Rank-Rank Scatter Plot

Outputs:

- `Auto_drug_reversal_rank_rank_scatter.pdf`
- `Auto_drug_reversal_rank_rank_scatter.png`
- raw point table: `Auto_drug_reversal_rank_rank_points.csv`

Generation:

The script inner-joins ASGARD and the secondary pipeline by:

- `state`
- `drug_key`

and plots:

- x-axis = ASGARD rank
- y-axis = secondary pipeline rank

Both axes are reversed so the best ranks sit near the lower-left corner of the panel.

The top-100 cutoffs are drawn as dashed lines at rank 100.

The script highlights:

- compounds in the top 100 of both methods
- top consensus compounds, with labels

Interpretation:

- lower-left compounds are strong by both methods
- compounds near one axis but far from the other are method-specific hits
- the density of points in the lower-left quadrant shows whether the two methods are converging on the same chemistry

## 13.3 Mechanism Target Dot Plot

Outputs:

- `Auto_drug_reversal_mechanism_target_dotplot.pdf`
- `Auto_drug_reversal_mechanism_target_dotplot.png`
- source data: `Auto_drug_reversal_top5_target_expression.csv`

Generation:

1. For each state, the script selects up to five top consensus compounds.
2. Preference is given to compounds that already have target annotations.
3. The target genes are parsed from:
   - `asgard_target_genes`
   - `scdrugprio_target_genes` or fallback target annotations
4. The script loads the five-state Seurat object.
5. For each target gene and each state, it computes:
   - mean log-normalized expression
   - fraction of cells expressing the gene
6. A dot plot is drawn where:
   - x-axis = expression state
   - y-axis = target gene
   - dot color = mean expression
   - dot size = fraction of cells expressing
   - an outline marks the intended target state

Interpretation:

- a target gene with a large, dark dot in the target state means the drug’s known target is actually expressed in that tumour state
- a target gene with uniformly absent or negligible expression suggests weaker mechanistic plausibility in the tumour cells
- strong target-state expression plus good consensus ranking is the main rationale for shortlisting a compound for validation

## 14. Output Inventory and How to Read Each File

### 14.1 Input and DEG files

- `Auto_drug_reversal_metadata.csv`
  - one row per retained cell
  - use to verify the exact cell universe

- `Auto_drug_reversal_state_sample_cell_counts.csv`
  - one row per sample-state combination
  - use to check per-state sample support

- `Auto_drug_reversal_degs_all_states.csv.gz`
  - full fresh or fallback state-vs-rest DEG table
  - use for any later re-ranking or sensitivity analysis

- `Auto_drug_reversal_signature_top150.csv`
  - compact exported up/down signatures
  - this is the best single file for inspecting the actual disease signatures submitted downstream

- `Auto_asgard_anchor_gene_diagnostic.csv`
  - genes with small effect size, broad prevalence, and non-significant state difference
  - use to inspect ASGARD “anchor”-like genes that are stable across malignant states

### 14.2 ASGARD outputs

- `Auto_asgard_raw_drug_identification.rds`
  - raw ASGARD return object
  - useful for debugging version-specific output structure

- `Auto_asgard_ranked_drugs.csv`
  - standardized ASGARD ranking
  - primary transcriptomic reversal table

- `Auto_asgard_status.csv`
  - whether the ASGARD run succeeded

### 14.3 scDrugPrio outputs

- `Auto_scdrugprio_ranked_drugs.csv`
  - standardized network-based ranking

- `Auto_scdrugprio_status.csv`
  - records whether the run used:
    - installed package or sourced functions
    - explicit resources or package defaults
    - gene-symbol or Entrez identifier space

This status file is important when interpreting biological credibility.

### 14.4 Local CMap fallback outputs

- `Auto_clue_ranked_drugs.csv`
  - direct local rank-matrix scoring output

- `Auto_clue_local_rankmatrix_status.csv`
  - confirms fallback scoring status

### 14.5 Consensus outputs

- `Auto_drug_reversal_top100_overlap_summary.csv`
  - best file for reporting overlap counts in slides

- `Auto_drug_reversal_consensus_drugs.csv`
  - final prioritized table combining ASGARD and the secondary pipeline

- `Auto_drug_reversal_rank_rank_points.csv`
  - raw data behind the rank-rank scatter plot

- `Auto_drug_reversal_top5_target_expression.csv`
  - raw data behind the mechanism target dot plot

## 15. Interpreting the Final Prioritized Compounds

The most credible candidates are compounds that satisfy as many of the following as possible:

1. High ASGARD rank.
2. High secondary-pipeline rank.
3. Present in the top-100 overlap.
4. Low `rank_sum` and low `mean_rank`.
5. Positive `state_specificity_margin`.
6. Known target genes detectable in the target state on the mechanism plot.

In practice:

- overlap alone is not enough
- a compound with overlap but no target expression is weaker mechanistically
- a compound with excellent target expression but no cross-method support is weaker computationally

The short list for wet-lab validation should therefore come from the intersection of:

- consensus support
- state specificity
- mechanistic target expression

## 16. Current Operational Caveats

### 16.1 ASGARD tissue proxy

The ASGARD reference currently uses the `stomach` LINCS tissue proxy because an oesophagus primary site was not available in the reference build path used here.

### 16.2 scDrugPrio resource quality matters

The scDrugPrio wrapper works technically with package example resources, but those example resources may not share the same perturbational universe as ASGARD/LINCS.

If top-100 overlap with ASGARD collapses to zero across all five states, the consensus script automatically replaces scDrugPrio with the direct local CMap fallback.

### 16.3 CLUE metadata endpoint

The current working metadata endpoint for target and MoA annotation is:

- `api/perts`

The legacy `rep_drug_target` and `rep_drug_moa` routes returned HTTP 500 during implementation and are not used in the current fallback annotator.

### 16.4 Fresh DEG mode is the production mode

The biologically correct production run is the fresh pooled Wilcoxon state-vs-rest DEG path:

- `AUTO_DRUG_DEG_MODE=findmarkers`
- `AUTO_FORCE_DRUG_DEGS=1` when refreshing the DEG cache

The `global` mode exists only so the wrapper stack can still run when a real DEG pass is temporarily blocked.

## 17. Practical Rerun Order

For a full accurate rebuild after any input or method change, the intended order is:

1. rerun `Auto_drug_reversal_inputs.R` with fresh DEGs
2. rerun `Auto_drug_reversal_asgard.R`
3. rerun `Auto_drug_reversal_scdrugprio.R`
4. rerun `Auto_drug_reversal_local_cmap.R`
5. rerun `Auto_drug_reversal_consensus_visuals.R`

This guarantees that:

- the DEG cache
- the top-150 signatures
- the ASGARD ranking
- the scDrugPrio ranking
- the local fallback ranking
- the final consensus figures

all come from the same underlying state definition and DEG derivation.

## 18. Current Refreshed Run (2026-04-23)

The current refreshed production-style run completed the input stage with:

- `AUTO_FORCE_DRUG_DEGS=1`
- `AUTO_DRUG_DEG_MODE=findmarkers`
- `AUTO_EXPORT_DRUG_MATRIX=0`

The resulting fresh DEG table contains:

- `94,387` state-vs-rest DEG rows

recorded in:

- `PDOs_outs/Auto_drug_reversal/Auto_drug_reversal_input_status.csv`

The refreshed final presentation outputs were rebuilt from:

- fresh `FindMarkers()` DEGs
- refreshed ASGARD ranking
- refreshed direct local CMap fallback ranking
- refreshed consensus visualization step forced to use `CLUE_FALLBACK_LOCAL` as the secondary ranking source

This was done because the scDrugPrio wrapper, although it completed successfully, still relies on the package example PPI/drug-target universe in the current environment and is therefore not the preferred secondary evidence layer for the current meeting-ready consensus figures.

Refreshed top-100 overlap counts in the current final consensus are:

- `Classic Proliferative`: `22`
- `Basal to Intest. Meta`: `30`
- `Stress-adaptive`: `21`
- `SMG-like Metaplasia`: `40`
- `3CA_EMT_and_Protein_maturation`: `32`

The current final consensus table contains:

- `145` overlapping drug-state rows

and the current final target-expression table contains:

- `290` rows across the five states

These values correspond to the current on-disk outputs in `PDOs_outs/Auto_drug_reversal/consensus/`.

## Predicted Reversion Visual Evidence

`analysis/cell_states/Auto_drug_reversal_predicted_reversion_visuals.R` adds presentation figures that explicitly test whether selected inhibitors look biologically plausible rather than merely appearing in ranked tables.

These plots are computational predictions. They do not use treated PDO RNA-seq data. The drug-side signal is inferred from the ASGARD/LINCS rank-matrix reference, scDrugPrio/DrugBank target annotations, and the scDrugPrio PPI network.

The waterfall plot ranks every drug within each method and state using rank evidence `-log10(rank / method universe)`. Final overlapping candidates are highlighted, so the audience can see whether the chosen inhibitors sit near the high-confidence end of the screen rather than being arbitrary labels.

The predicted anti-correlation scatter compares each gene's real PDO state-vs-rest logFC against the inferred opposing LINCS coordinate for the selected drug. Each dot is a state-defining gene. Genes with the largest absolute PDO logFC are highlighted. A strong negative trend in the opposing coordinate supports the hypothesis that the reference drug signature counteracts the malignant state signature.

The predicted state-flipping heatmap uses the same evidence in heatmap form. Rows are top state-defining genes. The first column is the real PDO malignant logFC. The second column is the predicted treatment opposition derived from the inverse centered LINCS rank coordinate. A visually convincing candidate should show a mirror-like pattern between the malignant column and predicted-treatment column.

The PPI targeted-hub overlay uses the scDrugPrio PPI and DrugBank action table. Nodes are disease-network genes colored by PDO logFC; larger nodes have higher PPI degree; diamond markers indicate annotated drug targets. This plot is intended to show whether an inhibitor targets connected/hub-like components of the malignant state network.

## scDrugPrio Directionality Safeguard

Network proximity is not sufficient for inhibitor nomination. A drug that activates an upregulated malignant-state target is a predicted mimic, not a predicted reversal, even if the target sits at the centre of the PPI network.

`Auto_drug_reversal_scdrugprio.R` therefore assigns every action-annotated drug-target pair to one of four classes:

- `counteracting`: an inhibitor/antagonist of an upregulated DEG, or an activator/agonist of a downregulated DEG.
- `mimicking`: an activator/agonist of an upregulated DEG, or an inhibitor/antagonist of a downregulated DEG.
- `mixed`: a multi-target drug with both counteracting and mimicking annotated DEG targets.
- `no_directional_target`: targets have no usable DEG direction in the current state.

The corrected scDrugPrio ranking excludes `mimicking` drugs by default with `AUTO_SCDRUGPRIO_EXCLUDE_MIMICS=1`, requires at least one counteracting DEG target by default with `AUTO_SCDRUGPRIO_REQUIRE_COUNTERACTION=1`, and ranks drugs by biological direction before PPI proximity. This prevents hub-proximal activators or no-direction target hits from being selected as top candidates solely because they strike central nodes.

`Auto_drug_reversal_scdrugprio_visuals.R` writes a refined network with explicit `target_direction`, `target_logfc`, and `target_status` fields. In that plot, a target labeled `No DEG information` should be treated as unproven directionally; it is not evidence of reversal.

For interpretability, the scDrugPrio visualization is allowed to show the pre-filtered network-proximity audit candidates. This is distinct from the final ranked scDrugPrio output: the audit plot is for mechanism review, while the final ranked table follows the pharmacological-action counteraction filter described above.
