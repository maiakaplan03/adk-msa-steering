# Phase-1 results

760 predictions: 152 edit recipes (7 columns × 3 mechanisms × 2 target directions ×
4 doses, minus the combinations excluded for the negative control) × 5 Boltz-2 seeds.

Each structure is summarised by `dTM = TM_Open − TM_Closed`. Positive means the
structure came out looking more like the open reference.

---

## 1. The seed moves the output 24× more than the edit does

Take the 760 numbers and average them two ways.

| grouping | n | variance of the means | what it measures |
|---|---:|---:|---|
| by edit recipe | 152 | **0.00040** | how much the MSA edit moves the output |
| by random seed | 5 | **0.00975** | how much the seed moves the output |

Ratio: **24.2×** for `dTM`, **39.8×** for the LID-Core angle.

`fig_seed_dominance.png` is the direct picture of this. Each violin contains 152
different edited alignments passed through a single seed. If the edits were doing the
work, each violin would be wide. Seeds 1, 2 and 3 are near-degenerate spikes — the same
structure emerges no matter which alignment goes in. Only seeds 4 and 5 show
edit-driven spread, which means the cell means are effectively being set by 2 of 5
draws.

This is the same failure mode Phase 0 flagged (within-cell std ≈ ensemble std,
ratio 0.93). Phase 1 inherits it with a tighter measurement.

## 2. The effects that exist are direction-blind

The hypothesis predicts opposite shifts for closed-targeted and open-targeted edits.
Observed at dose 100:

| effect | closed target | open target | opposite? |
|---|---:|---:|---|
| col 202, mechanism A, `dTM` | +0.066 | +0.064 | no |
| col 202, mechanism B, `dTM` | +0.068 | +0.068 | no |
| col 206, mechanism A, `dTM` | +0.011 | +0.017 | no |
| col 206, mechanism A, LID | 117° | 117° | no |
| col 64, mechanism C, `dTM` | +0.051 | **+0.062** | **yes, slightly** |

Aggregated over all 760 predictions, closed-targeted cells yield 31.8% open predictions
and open-targeted cells yield 32.4%. Whatever the edits are doing, the target direction
is not part of it.

The single exception — column 64, mechanism C, dose 100, open target — is the only place
the A/B/C × direction comparison points the way the hypothesis predicts, and the effect
(≈ +0.01 `dTM`) sits inside the cell standard error. It is a lead for a future
experiment, not a result.

## 3. A ceiling on how much steering was ever possible

Across all 760 predictions the NMP-Core angle spans **7.5°**. The same predictor on the
*unedited* MSA spans **37.5°**. The predictions fall into two tight clusters at roughly
(NMP 68°, LID 137°) and (NMP 67°, LID 102°); the references sit at (45°, 106°) for 1AKE
and (74°, 146°) for 4AKE. Neither cluster reaches either reference.

Under the absolute classification rule, Phase 1 gives 0% closed / 32.1% open / 67.9%
intermediate, against a null of 3.5% / 28.0% / 68.5%. Closed predictions disappear
entirely and the open fraction barely moves.

That narrowing is a property of the predictor on this protein, not of the edits — it
caps the headroom any input-side intervention could have used.

## 4. What the controls say

- **Negative-control column (50)**: null at all four doses in both directions, all
  within 0.001 of one another and within 0.02 of the null distribution. The pipeline
  is not manufacturing effects.
- **Mechanism B direction symmetry**: `|closed − open| < 0.01` for every (column, dose)
  except (136, 100), where the 0.031 gap traces to a single seed landing on the opposite
  attractor. With n = 5 that is sampling variance, not asymmetry.
- **Completeness**: all 760 structures present, all 214 CA atoms parsed, no dropouts.

## 5. Verdict on the hypothesis

The hypothesis had three parts. They resolve differently.

**"Conformational preference is positionally encoded in the MSA"** — partially supported.
The negative control is null while columns 136, 202 and 206 are not, so *where* you edit
matters. But the position-conditioned signal is 24× smaller than seed variance.

**"A small number of column edits can reproducibly steer the prediction"** — not
supported. The effects are direction-blind, and the state mix is unchanged between
closed-targeted and open-targeted edits.

**"Comparing mechanisms A/C against B distinguishes identity-reading from
co-evolution-reading"** — inconclusive. The underlying signal is too weak to support
the discrimination. Where the mechanisms do diverge, the pattern doesn't map onto either
side of the dichotomy: at column 136, A and B behave alike and C differs, which splits
by "preserves a distribution of residues" versus "installs one residue" rather than by
co-evolution versus identity.

One nuance worth carrying forward: at column 136 the active variable looks like
*anti-modal diversity* — the position holding many residues, none of them the modal one —
rather than co-evolution as originally framed.

## 6. What this says to do next

The design that produced this result also says clearly how to improve on it.

1. **More seeds, and report the median.** Seeds 1–3 return the same structure for any
   input; they are wasted reads. Going to ~20 seeds and switching to the median would
   roughly halve the cell-side standard error. Better still: identify which seeds are
   edit-responsive and weight accordingly.
2. **Multi-column edits.** Single-column edits move `dTM` by ≤ 0.07 here, while Phase 0's
   10+ column ensembles produced bimodal outcomes. Testing whether column effects compose
   is the natural next step and is cheap to add to `modifier_v2.py`.
3. **A second predictor.** Boltz's frozen NMP angle is a real ceiling. Running the same
   alignments through AlphaFold2 or another engine would separate "MSA editing doesn't
   work" from "this predictor is insensitive to MSA composition for this protein."
4. **Revisit column selection.** Columns 185 and 190 had the highest subfamily divergence
   on paper and produced the smallest effects; the columns that moved anything (136, 202,
   206) were not all top-divergence picks. The selection criterion did not pay off.
