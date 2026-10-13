# Evaluating the predictions

This page describes how the 1,296 ensemble predictions (main README, "Ensembles") were
scored and what they were compared against. The scoring was added during review; the
QA/QB map and the PCA of the pipeline itself are described in [`pipeline.md`](pipeline.md).

## Inputs

Three alignments from Guan et al. (2024), all sharing the same query (*E. coli* AdK) and
the same 214 match columns, with no lowercase insert columns:

| file | depth | what it is |
|---|---:|---|
| `MSA.a3m` | 2,651 | the full AdK sequence family; the base for every edit |
| `Closed500.a3m` | 501 | sequences with low frustration when threaded onto 1AKE (used only in the follow-up) |
| `Open100.a3m` | 101 | the same selection against 4AKE (used only in the follow-up) |

They are fetched, not copied into this repository — see `data/fetch_inputs.sh`.

## Prediction protocol

Every alignment was folded once with Boltz-2 v2.2.1, MSA subsampling on at 1,024 rows
(the default of the Python API; note the command-line flag caveat in
[`subsampling_audit.md`](subsampling_audit.md)). Because of that subsampling, folding
the same alignment twice can give different structures, so the reference for every
comparison is the **unedited MSA folded 314 times** under the same settings
(`results/sweep_evaluation/baseline_control_314.csv`).

## Metrics

**TM-score** to each crystal, over the 214 residue-aligned Cα atoms, with iterative
Kabsch superposition and the Zhang-Skolnick normalisation `d0 = 1.24·∛(L−15) − 1.8`,
L = 214. References are the Cα-only 1AKE and 4AKE structures of the upstream repository.

**Domain hinge angles**, between domain centroids (1-based residue numbers):

| domain | residues |
|---|---|
| NMP | 35–55 |
| CORE_anchor (vertex of the NMP angle) | 90–100 |
| hinge / corelid (vertex of the LID angle) | 115–124 |
| LID | 125–153 |
| CORE2 | 179–185 |
| CORE | everything else |

The implementation reproduces the values reported in the source paper: 1AKE at
NMP 45.1° / LID 106.4°, 4AKE at NMP 73.8° / LID 146.2°
(`evaluation/angle_analysis_plot_fixed.py`).

## How a prediction is classified

The leading rule is **absolute**: it asks whether a structure actually reached a state.

```
closed        if TM_Closed >= 0.85 and TM_Open   < 0.75
open          if TM_Open   >= 0.85 and TM_Closed < 0.75
both          if both >= 0.85                        (ambiguous)
intermediate  otherwise
```

A relative rule (`|TM_Closed − TM_Open| > 0.10`) is also computed in
`evaluation/sweep/05_verify_and_export.py`. It measures leaning rather than arrival and
gives higher closed/open fractions; it is reported alongside, never on its own. All
numbers in the README use the absolute rule.

## Scripts

`evaluation/sweep/` reads the Boltz-2 outputs of the four ensembles
(`ADK_RUNS`, default `../nuevos_MSA/fuente_MSA_original/results`) and writes
`results/sweep_evaluation/`:

| script | output |
|---|---|
| `01_param_sweep.py` | merged per-prediction table, state fractions per ensemble, parameter slices |
| `02_msa_diff.py` | which rows and columns each alignment changed (`col_edit_matrix.npy`) |
| `03_column_tendency.py` | per-column edit profile and its correlation with TM-score; per-domain summary |
| `04_plots.py` | `fig1`–`fig4` |
| `05_verify_and_export.py` | completeness checks; both classification rules |
| `06_steering_diagnostic.py` | within-cell vs between-cell spread |
| `07_target_distribution_audit.py` | checks on the edited columns' residue distributions |

`evaluation/mi_on_real_msa.py` tabulates the MI change of every ensemble alignment
(`results/sanity_check/mi_on_real_msa.csv`).

## What the alignment diffs show

- `break` leaves the amino-acid composition of every column unchanged — checked directly:
  at the most edited column (residue 101, 2,367 of 2,651 rows changed) the before/after
  residue counts are identical. Any effect of `break` must therefore go through the
  pairing of residues between columns.
- The ensemble that populates both ends (B−A `break`) is also the one that makes the
  fewest edits: on average 14,360 residue changes over 17.8 columns per alignment,
  against 22,859 over 43.1 for A−B `promote_joint` (`break` writes only one column of
  each pair, `promote_joint` writes both). Its edits are concentrated on CORE_anchor (residues 90–101, the vertex of the NMP angle), while the
  hinge is not touched at all.
- Every edited column is, by construction, a column of a state-specific contact, so these
  ensembles say nothing about columns outside the contact sets.
