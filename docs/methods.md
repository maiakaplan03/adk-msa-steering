# Methods

## Inputs

Three alignments, all sharing the same query (E. coli AdK) and the same 214 match
columns, with no lowercase insert columns:

| file | depth | what it is |
|---|---:|---|
| `MSA.a3m` | 2,651 | the full AdK sequence family; base for every edit |
| `Closed500.a3m` | 501 | sequences whose frustration is low when threaded onto 1AKE |
| `Open100.a3m` | 101 | same selection logic against 4AKE |

The two subfamilies come from Guan et al.'s frustration analysis and act as the target
distributions for directed edits. They are fetched, not vendored — see
`data/fetch_inputs.sh`.

## Contact pipeline (Phase 0)

Residue-residue contacts computed per conformation with a shadow map: heavy atoms only,
5.2 Å cutoff, 1.3 Å shadow radius (contacts occluded by an intervening atom are
discarded), sequence separation |Δseq| ≥ 3.

| | 1AKE (closed) | 4AKE (open) |
|---|---:|---:|
| contacts | 489 | 480 |
| conserved | 423 | 423 |
| conformation-specific | 66 | 57 |

**Known limitation.** 1AKE is ligand-bound (Ap5A) and the LID domain closes over that
ligand. Because the pipeline uses protein residues only, ligand-mediated LID closures
are invisible to it — the captured LID contacts (14-132, 15-132, 18-131) sit at the edge
of the interface. The pipeline therefore under-represents the LID's functional
importance.

## The editing recipe (Phase 1)

For each (column *j*, mechanism *T*, dose *d*, direction *D*):

```
R_j = { rows i : MSA[i][j] is not a gap, i != 0 }
S   = first floor(d · |R_j|) rows of R_j after a seed-0 shuffle

A — amplification
    P_target(j) = non-gap residue frequencies at column j in
                  Closed500.a3m  if D = closed
                  Open100.a3m    if D = open
    new residue := sample from P_target(j) with a per-row seed

B — break / anchor removal
    V = the multiset of residues currently at column j across R_j
    new residue := sample from V with a per-row seed
    (direction is inert here; both are run to confirm symmetry)

C — switch installation
    m_c = modal residue in Closed500 at column j
    m_o = modal residue in Open100   at column j
    if m_c == m_o:  target := 2nd-most-frequent residue of the target subfamily
    else:           target := m_c if D = closed else m_o
    new residue := target, identically for every row in S
```

Mechanism B changes no per-column amino-acid frequency by construction — it permutes
which row holds which residue. That is what makes the A/C-versus-B contrast a test of
composition against co-evolutionary pairing.

### Invariants, asserted at every call

- 214 columns exactly.
- Row 0 (the query) is never edited.
- Gaps are never introduced or removed; lowercase a3m inserts are treated as gaps.
- Row selection uses a fixed global seed of 0, so doses are nested:
  rows(25%) ⊂ rows(50%) ⊂ rows(75%) ⊂ rows(100%).
- Per-row RNG seeds are derived by SHA-256 from the row key, so the residue assigned to
  a given row is stable across dose tiers. A lower-dose alignment is therefore a literal
  subset of a higher-dose one at the level of individual substitutions, not merely of
  row selections.

Determinism is verifiable: `data/phase1_msa_checksums.sha256` holds the SHA-256 of all
152 alignments as used in the analysis.

## Column selection

| column | domain | why chosen | mechanism C viable? |
|---|---|---|---|
| 64 | CORE | high subfamily divergence; contact reorganization; V→I | yes, strong |
| 136 | LID | Phase-0 empirical sensitivity | weak (fallback target) |
| 185 | CORE2 | high divergence; E→S | yes |
| 190 | CORE | highest divergence overall (0.140); N→A | yes, strong |
| 202 | CORE | Phase-0 sensitivity; contact reorganization | weak (I↔L, conservative) |
| 206 | CORE | high divergence; F→R, aromatic to charged | yes, strong |
| **50** | NMP | **negative control** — lowest divergence, no contact change | skipped |

Columns 109, 155 and 172 were considered and dropped: their modal *and* second-modal
residues agree between subfamilies, so mechanism C has nothing to install.

## Metrics

**TM-score** to each reference, over 214 residue-aligned CA atoms, iterative Kabsch
superposition, Zhang-Skolnick normalization `d0 = 1.24·∛(L−15) − 1.8` with L = 214.
References are the CA-only 1AKE and 4AKE structures shipped with the upstream repository.
Implementation is numpy-only, in `src/score_boltz.py`.

**Domain hinge angles**, computed between domain centroids (1-based residue numbers):

| domain | residues |
|---|---|
| NMP | 35–55 |
| CORE_anchor (NMP angle vertex) | 90–100 |
| hinge / corelid (LID angle vertex) | 115–124 |
| LID | 125–153 |
| CORE2 | 179–185 |
| CORE | everything else |

The implementation reproduces the reference values reported in the source paper:
1AKE at NMP 45.1° / LID 106.4°, and 4AKE at NMP 73.8° / LID 146.2°.

**State classification** is documented in the README. The absolute rule
(0.85 / 0.75 thresholds) leads throughout; the relative rule (|dTM| > 0.10) is reported
alongside for Phase 0 only.

## Prediction protocol

All structures were predicted with Boltz-2 v2.2.1, MSA subsampling on at 1,024 rows
(the default behaviour of the Python API; note the CLI flag caveat in
`subsampling_audit.md`). Phase 1 used 5 seeds per edited alignment. The null
distribution is 314 repeats of the unedited alignment under the same settings.
