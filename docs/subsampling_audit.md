# Where the run-to-run variation actually comes from

Running the **identical** unedited alignment 314 times returns a wide spread of
structures — NMP-Core angle from 42° to 80°, LID-Core from 99° to 152°. The obvious
reading is that the model is expressing some conformational prior. It isn't.

## The audit

Checked against the Boltz-2 source (v2.2.1):

- **The only data-dependent stochastic operation in the trunk** is the MSA subsample,
  `src/boltz/model/modules/trunkv2.py`:
  `torch.randperm(msa.shape[1])[: self.num_subsampled_msa]`, gated on `subsample_msa`.
- **Eval-time dropout is inert.** `get_dropout_mask` in `model/layers/dropout.py`
  computes `dropout = dropout * training`; at inference `training = False`, so the rate
  is 0 and the mask is deterministically all ones. Dropout contributes nothing.
- **Featurization is deterministic.** `inferencev2.py` hardcodes
  `np.random.default_rng(42)`, so `--seed` does not change which sequences get featurized.
- **An upstream flag gotcha.** `--subsample_msa` is a click `is_flag` with no `default=`
  (`main.py:1022`), so it is **False when omitted from the command line** — even though
  the help text says "Default is True" and the Python `predict()` signature defaults to
  `True` (`main.py:1076`). Same input, two entry points, different behaviour.

So identical input producing different structures is caused by **which random 1,024 rows
get drawn from the alignment**, not by intrinsic model noise.

## Why this does not invalidate the null distribution

A control's job is to define a null: what spread of conformations comes back with no
editing at all, just from re-running the pipeline. That is exactly what 314 repeats of
the unedited alignment measures. The *mechanism* behind the spread does not change
*what the control measures* — subsampling noise is still the floor an editing effect has
to clear.

Two things follow from the data itself:

1. **The spread proves subsampling was on.** With subsampling off the trunk is
   deterministic and all 314 repeats would be identical. They are not, so the
   `randperm` path was active — the control is internally consistent with the finding.
2. **Depth rules out a hidden confound.** Subsampling only randomises when depth exceeds
   1,024. Every alignment here — unedited and all edited variants — has depth 2,651, so
   every condition draws 1,024 of 2,651. The edits change alignment *composition*, never
   depth-below-threshold, so no condition was accidentally deterministic.

## What it does change

**Vocabulary, not validity.** The ensemble from repeated same-input runs is not the
model expressing a conformational prior; it is the sensitivity of a single-structure
predictor to alignment composition, surfaced by random subsampling. That is the
AlphaFold-era MSA-subsampling technique for eliciting alternative conformations
(Del Alamo et al. 2022; AF-Cluster) — the right lineage to cite.

It also puts the editing experiment squarely in the right conceptual frame: if
conformation is driven by alignment composition, then deliberately changing composition
is the correct lever to pull. The framing is reinforced, not undermined.

## Caveats kept open

- **Provenance.** The unedited control's spread proves subsampling was on there; the
  edited batches should be confirmed the same way (`subsample_msa: true`,
  `num_subsampled_msa: 1024` in each run's `hparams.yaml`) so the comparison is provably
  matched rather than presumed.
- **Composition axis.** A direction in the pair representation that separates open from
  closed could partly be an alignment-composition axis that correlates with
  conformation. For the editing experiment that is expected; any claim that the
  direction is *purely* conformational needs a control checking whether the axis also
  separates runs by subsample summary statistics.
- **Reproducibility going forward.** Batches launched without an explicit `--seed` were
  entropy-seeded: diverse, but not reproducible. Pass a seed per run.
  `evaluation/check_seed_determinism.sh` distinguishes the three candidate sources of
  run-to-run variation (unfixed seed, RNG advancing within one process, and genuine
  GPU/bf16 nondeterminism across processes).
