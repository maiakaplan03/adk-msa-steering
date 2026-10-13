# Review notes

What was checked, changed or questioned when the pipeline was reviewed with Claude
(Anthropic). Each note says which step and which case it applies to, what was found, why
it matters, and what was done. Nothing here changes the alignments that were folded:
with its defaults, the code in `pipeline/` regenerates all 1,296 ensemble alignments byte
for byte.

## What was verified

- **The editor is reproducible.** `pipeline/04_generate_sweep.py --ensemble all`
  regenerates the 1,296 alignments of the four ensembles exactly
  (`data/sweep_msa_checksums.sha256`, 1,296/1,296 matching).
- **The contacts are reproducible.** `pipeline/01_contacts_shadow.py` on the upstream
  `1ake_A.pdb` gives the same 489 contacts as `results/contacts/`, and
  `pipeline/02_delta_contacts.py` reproduces the lost / gained / conserved / per-residue
  tables exactly.
- **The toy MI check is reproducible.** `pipeline/03_sanity_check_mi.py` gives the same
  picture as the presentation (`promote`, `promote_joint` 100% up; `break` 50–73% down).
- **The query row is never edited, no gap is added or removed**, in any of the 1,296
  alignments (`results/sweep_evaluation/msa_diff_summary.csv`).

## 1 · Shared columns are overwritten in order (step 3, `apply_contacts_in_msa`)

**Found.** All selected pairs edit the same column lists, one pair after another. If
residue 101 is in three selected pairs, its column is edited three times and keeps the
last edit; the "MI before" of the second and third pair is measured on a column that
was already changed. `degree` mode picks exactly the high-degree residues, so this
happens most where the edit is strongest. Example (A−B, frac 0.9, strength 0.9, seed 11,
min_seq_sep 8): 54 pairs are chosen and 26 residues appear in more than one of them.

**Why it matters.** The result depends on the order the pairs were drawn in, and an
earlier pair's edit can be partly undone by a later one.

**Done.** The original behaviour is kept as the default (`overlap="sequential"`),
because it is what was folded. `overlap="lock"` (`--overlap lock` in the sweep script)
lets each column be written by only one pair; later pairs that would overwrite it are
skipped and marked `skipped=True` in the report. In the example above, `lock` skips 25
of 54 pairs for `promote_joint` (which writes both columns) and 10 of 54 for `break`
(which writes only column B).

## 2 · `promote_joint` lowers MI at medium and high strength (steps 3–4, real alignment)

**Found.** On the toy alignment `promote_joint` always raised MI. On the real AdK
alignment, over the 648 `promote_joint` alignments of the ensembles:

| strength | share of alignments where mean ΔMI > 0 | mean ΔMI (bits) |
|---|---:|---:|
| 0.3 | 81–85% | +0.003 to +0.004 |
| 0.6 | 22–58% | −0.004 to +0.002 |
| 0.9 | 4–6% | −0.013 to −0.022 |

(ranges are over the A−B and B−A ensembles; `break` lowered MI in 100% of its 648.)

**Why.** `promote_joint` draws new pairs from the joint distribution after raising it
to `gamma` (2–4). That concentrates almost all the probability on the most common pair,
so at high strength most edited rows get the same two residues. In a typical
strength-0.9 alignment (A−B, frac 0.9, gamma 3) the edited columns fall from 1.45 to
0.60 bits of entropy and their most common residue goes from 65% to 85% of rows. MI
can never exceed the entropy of either column, so a column that is nearly constant
carries almost no MI. The toy alignment did not show this because its coupled pair was
built from five different residue pairs spread evenly.

**Why it matters.** At medium and high strength, the `promote_joint` ensembles test
"make these columns conserved" rather than "make these columns coevolve". This is the
most likely reason they stayed close to the null.

**Suggestion.** For the next sweep, keep `promote_joint` at low strength, or use
`gamma` = 1, and check ΔMI on the real alignment before folding
(`evaluation/mi_on_real_msa.py`).

## 3 · `gamma`, `alpha` and the MI itself (step 3)

- `gamma` (`power_gamma`) is not a detail: it decides how far the sampling distribution
  is pushed towards its most common outcome, and it was swept (2, 3, 4). It should be
  named as a parameter of the method.
- `alpha` behaves slightly differently in the two promote actions: in `promote` it is
  added to every residue seen in column B for each residue of A; in `promote_joint` it is
  added only to pairs that already occur, so a pair never seen in the alignment can
  never be created.
- MI is raw: no correction for finite sample size, sequence redundancy or the
  "average product correction" (APC) used in coevolution papers. That is fine for its
  job here — comparing the same pair before and after an edit — but the values are not
  comparable with published MI or DCA scores.

## 4 · Which predictor (step 5)

Everything was folded with **Boltz-2** (v2.2.1). Labels saying "AF2" in early code and
tables were leftovers from the original plan and were corrected (for example the column
`af2_seed` → `boltz_seed`).

## 5 · Why the same MSA gives different structures (step 5)

Boltz-2 randomly picks 1,024 of the 2,651 sequences of the alignment on every run, so
folding an unedited alignment repeatedly gives a spread of structures (NMP angle 42–80°,
LID angle 99–152° over 314 runs). This was confirmed in the Boltz-2 source code
([`subsampling_audit.md`](subsampling_audit.md)). Two consequences:

- any edited ensemble has to be compared with **the unedited alignment folded many
  times**, not with one unedited prediction — hence the 314-run null;
- every alignment here has the same depth (2,651), so every condition is subsampled the
  same way and the comparison is fair.

Still to confirm: that the edited batches also ran with `subsample_msa: true` and
`num_subsampled_msa: 1024` (in each run's `hparams.yaml`).

## 6 · Scoring the ensembles (step 5)

The presentation evaluated structures with QA/QB and PCA. For the ensembles, TM-score to
each crystal and the two hinge angles were added, with an **absolute** rule (a
prediction is closed or open only if it reaches TM ≥ 0.85 to that crystal and < 0.75 to
the other). A relative rule (`|ΔTM| > 0.10`) is reported alongside only, because it
calls "closed" anything that leans closed, even if it is far from 1AKE
([`evaluation.md`](evaluation.md)).

## 7 · One prediction per alignment (step 5, ensembles)

Each of the 1,296 alignments was folded once. The state fractions in the README are
therefore over 324 different alignments, not over repeated runs of the same one, and
part of each fraction is subsampling noise (note 5). The clearest pattern — `break`
moving predictions away from the state whose contacts were broken — is large compared
with that noise (closed 3.5% → 14.2%, open 28.0% → 35.8%), but it is one run per
alignment and one protein.

## 8 · The Ap5A ligand (steps 0–2)

1AKE is solved with the inhibitor Ap5A, and the LID domain closes over it. Contacts are
computed on protein atoms only, so contacts that exist through the ligand are
invisible: the LID contacts that are captured (14–132, 15–132, 18–131) sit at the edge
of the interface. The Δcontact sets therefore under-represent the LID.

## 9 · The noise window in the PCA (step 5.2)

The notebook comment said "ruido: −100 ≤ PC1 ≤ −40", but the code marked
`40 ≤ PC1 ≤ 150` — the cluster on the right of the first PCA plot, which was the one
removed (171 structures). The window is now a parameter (`--noise-pc1-min`,
`--noise-pc1-max`, default 40 and 150). The sign and position of PC1 change from one
data set to another, so the window has to be read off the first plot every time.

## 10 · Slide vs. code parameters (step 4)

The presentation lists frac 0.2–0.8 and strength 0.4–0.8 with `promote` and
`promote_joint`; the ensembles that were scored used frac and strength 0.3/0.6/0.9 with
`promote_joint` and `break`. Both grids are now named presets of the sweep script
(`presentation`, `ensembles`) so it is clear which numbers belong to which experiment.

## 11 · Ties in the top-10 lists (step 2)

`delta_top_gainers_with_partners.csv` and `delta_top_losers_with_partners.csv` take the
ten residues with the largest degree change. Several residues tie (for example many at
+2), and the order among ties depends on the pandas version, so a re-run can list a
different tenth residue. The full per-residue table is identical.
