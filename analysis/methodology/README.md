# PDO Analysis Methodology Index

Methodology files are required for active workflows with non-trivial scientific
choices, thresholds, statistical models, or external references. Simple direct
exports and plot-only scripts use `Methodology: none` in their registry header.

## Active methodology

- `metaprograms/centred/Auto_centred_metaprogram_refinement_methodology.md`:
  centred rank selection, parent QC, splitting/merging, and final MP QC.
- `metaprograms/centred/Auto_centred_ordering_and_state_methodology.md`:
  final MP order, UCell, state thresholds, 3CA comparison, and TCGA survival.
- `cell_states/Auto_centred_five_state_marker_methodology.md`: current
  centred-state marker discovery and recurrence thresholds.
- `cell_states/Auto_five_state_surface_marker_methodology.md`: surface
  annotation, topology rules, margins, and FACS ranking.
- `cell_states/centred_refined_pdo_flot_matched_response_methodology.md`:
  matched-FLOT composition, MP, pathway, and paired edgeR analyses.
- `cell_states/Auto_PDO_treatment_scenic_methodology.md`: matched treatment
  SCENIC filtering, AUC/RSS comparisons, and patient-pair interpretation.
- `cell_states/Auto_PDO_scRef_centred_state_concordance_methodology.md`:
  independent PDO/scRef centred-state concordance.
- `metaprograms/Auto_3CA_state_correlation_crossdata_methodology.md`,
  `Auto_3CA_pseudobulk_correlation_crossdata_methodology.md`, and
  `Auto_parse_external_centred_mp_timecourse_methodology.md`: current
  cross-dataset comparisons.
- `cnv/`, `cnv/wes_subclone/`, `demultiplex/`, and
  `demultiplex/strelka/`: optional diagnostic/tool workflows.

## Legacy methodology

Files prefixed `legacy_`, plus the drug-reversal, trajectory, old clinical,
old enrichment, old uncentred metaprogram, and old state-workflow methodology
documents, describe provenance only. They must not be interpreted as the
current centred route.

Every active methodology must state the aim, exact input/reference sources,
algorithmic steps, consequential thresholds and their rationale, statistical
unit and limitations, persistent outputs, cache/replot behavior, and downstream
dependencies.

