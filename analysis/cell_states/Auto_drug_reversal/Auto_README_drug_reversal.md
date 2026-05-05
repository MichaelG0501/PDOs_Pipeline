# Auto Drug Reversal Workflow

This folder is the canonical, reproducible layout for the PDO malignant-state drug-reversal analysis. It consolidates the active input preparation, ASGARD, scDrugPrio, CLUE/CMap fallback, local CMap, reference-download, and visualization scripts.

The original root-level and `analysis/cell_states/` copies were left in place for file safety. New work should use the copies in this folder.

## Active Workflow Files

| Purpose | R script | PBS / shell wrapper |
| :--- | :--- | :--- |
| Fresh per-state DEG/signature inputs | `Auto_drug_reversal_inputs.R` | `Auto_run_drug_reversal_inputs.sh` |
| ASGARD reversal scoring | `Auto_drug_reversal_asgard.R` | `Auto_run_drug_reversal_asgard.sh` |
| scDrugPrio transcriptomic/network prioritization | `Auto_drug_reversal_scdrugprio.R` | `Auto_run_drug_reversal_scdrugprio.sh` |
| CLUE API fallback | `Auto_drug_reversal_clue_fallback.R` | `Auto_run_drug_reversal_clue_fallback.sh` |
| Local CMap/L1000 fallback | `Auto_drug_reversal_local_cmap.R` | `Auto_run_drug_reversal_local_cmap.sh` |
| ASGARD reference preparation | `Auto_prepare_asgard_reference.R` | `Auto_build_asgard_reference.sh` |
| Combined method visualizations | `Auto_drug_reversal_method_visuals.R` | `Auto_run_drug_reversal_method_visuals.sh` |
| scDrugPrio network visualizations | `Auto_drug_reversal_scdrugprio_visuals.R` | `Auto_scdrugprio_viz.sh` |
| Conda environment | `Auto_drug_reversal_environment.yml` | `Auto_setup_drug_reversal_env.sh` |
| ASGARD L1000 reference download | - | `Auto_download_asgard_l1000_reference.sh` |

## Execution Order

Run from the project root:

```bash
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD
```

1. Prepare or update the software environment if needed:

```bash
qsub analysis/cell_states/Auto_drug_reversal/Auto_setup_drug_reversal_env.sh
```

2. Download and prepare reference resources if missing:

```bash
qsub analysis/cell_states/Auto_drug_reversal/Auto_download_asgard_l1000_reference.sh
qsub analysis/cell_states/Auto_drug_reversal/Auto_build_asgard_reference.sh
```

3. Regenerate fresh per-state DEG signatures:

```bash
qsub analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_inputs.sh
```

4. Run the three reversal methods:

```bash
qsub analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_asgard.sh
qsub analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_scdrugprio.sh
qsub analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_local_cmap.sh
```

5. Run CLUE API fallback only when required:

```bash
qsub analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_clue_fallback.sh
```

6. Generate presentation visualizations:

```bash
qsub analysis/cell_states/Auto_drug_reversal/Auto_run_drug_reversal_method_visuals.sh
qsub analysis/cell_states/Auto_drug_reversal/Auto_scdrugprio_viz.sh
```

## External Reference Locations

Large reference files are stored outside the repo:

- scDrugPrio PPI: `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/ppi.txt`
- Drug-target table: `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/all_drug_targets_drug_bank.txt`
- ASGARD/L1000 references: `/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000/`

The wrappers set or consume the relevant environment variables. If a resource is missing, update the reference folder rather than writing large databases into the project repository.

## Main Outputs

Outputs are written under:

```text
PDOs_outs/Auto_drug_reversal/
```

Expected suboutputs include the fresh DEG signatures, ASGARD ranked drugs, scDrugPrio ranked/audit tables, CLUE/local CMap ranked tables, and method-specific visualization PDFs/PNGs.

## State Order

Any plot or table touching finalized malignant state order should use:

1. Classic Proliferative
2. Basal to Intest. Meta
3. SMG-like Metaplasia
4. Stress-adaptive
5. Immune Infiltrating, only when present in scATLAS-style state sets
6. 3CA_EMT_and_Protein_maturation

Metaprogram order should follow this state order, and within each state should use `mp_tree_order` when available.
