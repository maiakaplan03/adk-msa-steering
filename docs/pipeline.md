# The pipeline, step by step

Maia Kaplan's pipeline, in the order of the TP presentation. Every step is a script in
`pipeline/`; the original notebook is `notebooks/TP_final_demo.ipynb`. Where the review
changed or questioned something, the note number in [`review_notes.md`](review_notes.md)
is given in brackets.

## 0 · Crystals — `00_get_crystals.py`

Downloads 1AKE and 4AKE from the PDB and keeps chain A, standard amino acids only.
Waters, ions and the Ap5A ligand of 1AKE are removed [note 8].

## 1 · Contacts — `01_contacts_shadow.py`

Residue–residue contacts with the **Shadow** criterion (Noel et al. 2012, the contact map
used by structure-based models):

- two heavy atoms are neighbours if they are within **C = 5.2 Å**;
- the pair is discarded if a third atom lies within **S = 1.3 Å** of the segment joining
  them (it "shadows" the contact); for pairs exactly at the minimum separation the shadow
  radius is reduced to 0.5 Å;
- residues closer than **3** positions in sequence are excluded;
- a residue pair is in contact if any of its atom pairs is, and the shortest such
  distance is kept.

Result: 489 contacts in 1AKE and 480 in 4AKE (`results/contacts/*_residue_contacts_shadow.csv`
and the 214 × 214 matrices). Re-running the script on the upstream `1ake_A.pdb` gives the
same 489 contacts.

## 2 · Δcontacts — `02_delta_contacts.py`

With A = 1AKE (closed) and B = 4AKE (open):

| set | file | n | meaning |
|---|---|---:|---|
| lost, A − B | `delta_contacts_lost_A_minus_B.csv` | 66 | closed-specific contacts |
| gained, B − A | `delta_contacts_gained_B_minus_A.csv` | 57 | open-specific contacts |
| conserved | `delta_contacts_conserved.csv` | 423 | in both, with the distance change |

`delta_contacts_per_residue.csv` lists, for every residue, its degree in A and B and its
gained and lost partners; `delta_top_*` holds the ten residues that gain or lose most
[note 11].

## 3 · The MSA editor — `msa_editor.py`

**Input:** the original MSA and one Δcontact list. **Output:** a synthetic MSA with the
same rows, columns and gaps, and a per-pair report of MI before and after.

Lowercase a3m insert characters are turned into gaps, residue numbers of the query
(row 0) are mapped to MSA columns, and the query row is never edited.

### Choosing pairs (mode)

From the contact list, pairs closer than `min_seq_sep` in sequence are removed, then each
remaining pair (i, j) gets a weight:

| mode | weight | idea |
|---|---|---|
| `uniform` | 1 | every pair equal |
| `degree` | deg(i) + deg(j) | residues in many contacts are more "connected" |
| `dist` | 1 / d(i, j) | closer contacts count more |
| `hybrid` | (deg(i) + deg(j)) / d(i, j) | connected and close |
| `random` | — | null: random residue pairs at the same separation, as many as in the list |

`k = round(frac × number of pairs)` pairs are drawn without replacement with probability
proportional to the weight. Degrees are counted inside the filtered contact list.

### Editing each pair (action)

For each chosen pair the editor takes the rows where neither column has a gap (minus the
query) and picks `round(strength × that number)` of them at random.

- **`promote`** — for each picked row with residue *a* in column A, column B is replaced
  by a residue drawn from P(B | A = a), estimated from the whole alignment with
  pseudocount `alpha` over the residues seen in B, then raised to the power `gamma` and
  renormalised. B is pushed to "agree" with A.
- **`promote_joint`** — the observed residue pairs (a, b) are counted, `alpha` is added,
  the distribution is raised to `gamma` and renormalised, and each picked row receives a
  new pair (a, b) drawn from it in **both** columns [note 2].
- **`break`** — the residues of column B in the picked rows are shuffled among those
  rows. Column B keeps exactly the same amino-acid composition; only which residue sits
  next to which residue of A changes. This is the negative control.

`gamma` > 1 sharpens the distribution towards its most common outcome; `gamma` = 1 leaves
it as estimated [note 3].

Pairs are processed in the order they were drawn, on shared column lists [note 1].

### 3.3 · Coevolution metric: mutual information

For two columns A, B, over the rows with no gap in either:

```
MI(A, B) = Σ_a,b p(a, b) · log2( p(a, b) / (p(a) · p(b)) )
```

Higher MI means the residues of the two columns depend more on each other. It is
recorded for every edited pair before and after the edit [note 3].

### 3.4 · Control with a synthetic query and MSA — `03_sanity_check_mi.py`

A synthetic alignment of 120 sequences × 40 positions where position 5 determines
position 30 in 70% of the rows (K→E, R→D, A→A, V→V, G→G), plus 13 contacts with a hub on
residue 5. Each action × mode (uniform, degree, hybrid) is run with 40 seeds at
frac 0.2, strength 0.6, alpha 0.05, gamma 4:

| action | expected | seeds in the expected direction |
|---|---|---|
| `promote` | MI up | 100% in all modes |
| `promote_joint` | MI up | 100% in all modes |
| `break` | MI down | 50% (uniform) to 70–73% (degree, hybrid) |

`break` is weaker because shuffling a subset of rows can, by chance, leave most pairs
intact; on the real AdK alignment it lowered MI in 100% of alignments
(`results/sanity_check/mi_on_real_msa.csv`).

## 4 · The sweep — `04_generate_sweep.py`

Every combination of the grid is generated and copied into one folder, with the
parameters in the file name:

```
msa_degree_frac0.30_str0.30_seed11__ACT-promote_joint__MODE-degree__A-0_05__G-2__MINSEP-12.a3m
```

| preset | actions | frac | strength | seeds | alpha | gamma | min_seq_sep | n |
|---|---|---|---|---|---|---|---|---:|
| `presentation` | promote, promote_joint | 0.2, 0.4, 0.6, 0.8 | 0.4, 0.6, 0.8 | 11, 22, 33 | 0.05, 0.1 | 2, 3, 4 | 8, 12 | 864 |
| `ensembles` | promote_joint *or* break | 0.3, 0.6, 0.9 | 0.3, 0.6, 0.9 | 11, 22, 33 | 0.05, 0.1 | 2, 3, 4 | 8, 12 | 324 per ensemble |

All in `degree` mode. The four ensembles are A−B and B−A, each with `promote_joint` and
with `break` (the folders called `control_negativo`). The presentation-era file names
also carried `__DF-0` (distance filter off) [note 10].

## 5 · Folding with Boltz-2

Each synthetic MSA is given to Boltz-2 with the AdK query, producing one structure per
alignment. See [`evaluation.md`](evaluation.md) for the settings and
[`subsampling_audit.md`](subsampling_audit.md) for why the same alignment can give
different structures [notes 5–6].

## 5.1 · QA vs QB — `05a_qa_qb.py`

The Shadow contacts of every predicted structure are computed with the same parameters
as step 1. For each prediction:

- **QA** = how many of the 66 closed-specific contacts it has;
- **QB** = how many of the 57 open-specific contacts it has.

A perfect closed structure sits at (66, 0) and a perfect open one at (0, 57); a path
between the states would run along the diagonal between them.

## 5.2 · PCA — `05b_pca.py`, `05c_pca_paths.py`

1. Every structure (predictions plus the two crystals) is superposed on 4AKE using its
   Cα atoms (Kabsch), and the 214 × 3 Cα coordinates are flattened into one row.
2. PCA of those rows; PC1 and PC2 are plotted with the crystals marked.
3. Structures in a separate cluster along PC1 are listed as "ruido" (noise) and the PCA
   is repeated without them [note 9].
4. A k-nearest-neighbour graph (k = 18) is built in the PC1–PC2 plane, its minimum
   spanning tree is computed, and the tree path is followed from the prediction closest
   to 1AKE to the prediction closest to 4AKE (and to a second reference).
5. Each path is thinned to about 15 structures spread evenly along its length, without
   the crystals — the candidate intermediates.

## Next steps (from the presentation)

- Build the full contact universe, to relate the formation of each individual contact to
  the structure Boltz-2 produces.
- Apply the algorithm to proteins with only one of the two conformations crystallised.
