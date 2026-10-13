"""
msa_editor.py — Maia Kaplan's pair-based MSA editor (TP final, "BLOQUE 1 — CORE FUNCTIONS").

Edits PAIRS of MSA columns that correspond to residue contacts which differ between
two conformations (AdK: 1AKE closed vs 4AKE open), and measures the mutual information
(MI) of each pair before and after the edit.

Actions (her presentation, section 3.1):
  promote        rows i of columns A,B: resample B from P(B | A=a)            -> MI up
  promote_joint  rows i of columns A,B: resample (A,B) from the joint P(A,B)   -> MI up
  break          rows i of column B: shuffle the residues among the rows      -> MI down
                 (negative control: column composition unchanged, pairing destroyed)

Pair-selection modes (section 3.2): uniform | degree | dist | hybrid  (+ random as a null)
Knobs: frac (fraction of pairs edited), strength (fraction of valid rows edited per pair),
seed, alpha (pseudocount), gamma (power_gamma: sharpening exponent applied to the
conditional / joint distribution before sampling).

This file is the notebook code moved into a module, unchanged in behaviour, with two
additions made during review (both OFF by default, so the defaults regenerate the
original 1,296 sweep alignments byte-for-byte — see data/sweep_msa_checksums.sha256):

  * overlap="lock" in apply_contacts_in_msa(): by default ("sequential") every selected
    pair edits the same shared column lists in order, so a column that belongs to
    several selected pairs keeps whichever edit came last, and the MI "before" of a later
    pair is measured on already-edited columns. With "lock", a column can be written by
    at most one pair; later pairs that would overwrite it are skipped and flagged in the
    report. See docs/review_notes.md, note 1.
  * write_reports=True in make_progressive_msas_weighted() writes the per-pair MI report
    (report_<tag>.csv) next to each alignment. Reports never change the .a3m output.

The original Spanish comments are kept.
"""
import os, math, json, re, shutil, random
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from typing import List, Optional, Tuple

# ---------- A3M / gaps ----------
GAP_CHARS = set(["-", "."])
LOWER_RE = re.compile(r"[a-z]")

def is_gap_char(ch: str) -> bool:
    return (ch in GAP_CHARS) or (len(ch) == 1 and ch.islower())

def read_fasta_or_a3m(path: str):
    headers, seqs = [], []
    with open(path) as f:
        h, s = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    seqs.append("".join(s))
                h = line[1:]
                headers.append(h)
                s = []
            else:
                s.append(line)
        if h is not None:
            seqs.append("".join(s))
    return headers, seqs

def write_fasta(path: str, headers: List[str], seqs: List[str]):
    with open(path, "w") as f:
        for h, s in zip(headers, seqs):
            f.write(f">{h}\n{s}\n")

def normalize_a3m(headers, seqs):
    # minúsculas (inserts A3M) -> '-' para no contaminar MI/acciones
    out = [LOWER_RE.sub("-", s) for s in seqs]
    return headers, out

# ---------- mapping ref (nogap residue index -> MSA column) ----------
def build_ref_nogap_to_col(ref: str, one_based=True):
    mapping = {}
    res = 0
    for col, ch in enumerate(ref):
        if not is_gap_char(ch):
            mapping[(res + 1) if one_based else res] = col
            res += 1
    return mapping

# ---------- MI ----------
def column_MI(colA, colB, ignore_gaps=True) -> float:
    pairs = []
    for a, b in zip(colA, colB):
        if ignore_gaps and (is_gap_char(a) or is_gap_char(b)):
            continue
        pairs.append((a, b))
    if not pairs:
        return 0.0
    n = len(pairs)
    cA = Counter(a for a, _ in pairs)
    cB = Counter(b for _, b in pairs)
    cAB = Counter(pairs)
    mi = 0.0
    for (a, b), nab in cAB.items():
        pa = cA[a] / n
        pb = cB[b] / n
        pab = nab / n
        mi += pab * math.log(((pab / (pa * pb + 1e-12)) + 1e-12), 2)
    return mi

def _indices_non_gap_both(colA, colB):
    return [i for i, (a, b) in enumerate(zip(colA, colB)) if (not is_gap_char(a) and not is_gap_char(b))]

# ---------- P(B|A) ----------
def build_conditional_P_B_given_A(colA, colB, alpha=0.1, power_gamma=3.0):
    pair_counts = defaultdict(Counter)
    A_support = set()
    B_support = set()
    for a, b in zip(colA, colB):
        if is_gap_char(a) or is_gap_char(b):
            continue
        pair_counts[a][b] += 1
        A_support.add(a)
        B_support.add(b)
    if not A_support or not B_support:
        return {}

    B_list = sorted(B_support)
    cond = {}
    for a in A_support:
        raw = np.array([pair_counts[a][b] + alpha for b in B_list], dtype=float)
        raw /= raw.sum()
        if power_gamma is not None and power_gamma != 1.0:
            raw = raw ** power_gamma
            raw /= raw.sum()
        cond[a] = (B_list, raw)
    return cond

# ---------- actions ----------
def promote_pair_columns(colA, colB, strength=0.6, rng_py=None, rng_np=None,
                         alpha=0.1, power_gamma=3.0, ref_index=None):
    if rng_py is None:
        rng_py = random.Random()
    if rng_np is None:
        rng_np = np.random.default_rng()

    idxs = _indices_non_gap_both(colA, colB)
    if ref_index is not None and 0 <= ref_index < len(colA):
        idxs = [i for i in idxs if i != ref_index]
    if not idxs:
        return

    cond = build_conditional_P_B_given_A(colA, colB, alpha=alpha, power_gamma=power_gamma)
    if not cond:
        return

    k = int(round(strength * len(idxs)))
    if k <= 0:
        return
    sel = rng_py.sample(idxs, min(k, len(idxs)))

    for i in sel:
        a = colA[i]
        if a in cond:
            B_vals, B_probs = cond[a]
            j = rng_np.choice(len(B_vals), p=B_probs)
            colB[i] = B_vals[j]

def break_pair_columns(colA, colB, strength=0.6, rng_py=None, ref_index=None):
    if rng_py is None:
        rng_py = random.Random()
    idxs = _indices_non_gap_both(colA, colB)
    if ref_index is not None:
        idxs = [i for i in idxs if i != ref_index]
    if len(idxs) <= 1:
        return

    k = int(round(strength * len(idxs)))
    if k <= 1:
        return

    sel = rng_py.sample(idxs, min(k, len(idxs)))
    buf = [colB[i] for i in sel]
    rng_py.shuffle(buf)
    for i, v in zip(sel, buf):
        colB[i] = v

# def promote_by_joint(colA, colB, strength=0.6, rng_np=None, ref_index=None):
#     if rng_np is None:
#         rng_np = np.random.default_rng()

#     idxs = _indices_non_gap_both(colA, colB)
#     if ref_index is not None:
#         idxs = [i for i in idxs if i != ref_index]
#     if len(idxs) == 0:
#         return

#     pairs = [(colA[i], colB[i]) for i in idxs]
#     vals, cnts = np.unique(pairs, axis=0, return_counts=True)
#     probs = cnts / cnts.sum()

#     k = int(round(strength * len(idxs)))
#     if k <= 0:
#         return

#     sel_positions = rng_np.choice(len(idxs), size=min(k, len(idxs)), replace=False)
#     sampled = vals[rng_np.choice(len(vals), size=len(sel_positions), p=probs)]
#     for pos, (a_new, b_new) in zip(sel_positions, sampled):
#         colA[idxs[pos]] = a_new
#         colB[idxs[pos]] = b_new

def promote_by_joint(colA, colB, strength=0.6, rng_np=None, ref_index=None,
                     alpha=0.1, power_gamma=3.0): # Add alpha/gamma here
    if rng_np is None:
        rng_np = np.random.default_rng()

    idxs = _indices_non_gap_both(colA, colB)
    if ref_index is not None:
        idxs = [i for i in idxs if i != ref_index]
    if len(idxs) == 0:
        return

    # 1. Get unique pairs and counts
    pairs = [(colA[i], colB[i]) for i in idxs]
    vals, cnts = np.unique(pairs, axis=0, return_counts=True)

    # 2. Apply Laplace Smoothing (Alpha)
    # This prevents zero-probability issues
    probs = cnts.astype(float) + alpha
    probs /= probs.sum()

    # 3. Apply Sharpening (Gamma)
    if power_gamma is not None and power_gamma != 1.0:
        probs = probs ** power_gamma
        probs /= probs.sum()

    k = int(round(strength * len(idxs)))
    if k <= 0:
        return

    sel_positions = rng_np.choice(len(idxs), size=min(k, len(idxs)), replace=False)
    sampled = vals[rng_np.choice(len(vals), size=len(sel_positions), p=probs)]

    for pos, (a_new, b_new) in zip(sel_positions, sampled):
        colA[idxs[pos]] = a_new
        colB[idxs[pos]] = b_new

# ---------- contact weights + filters ----------
DIST_COLUMN_PRIORITY = ["dist_B", "dist", "Dist", "DIST", "dist_A"]

def _pick_distance_column(df: pd.DataFrame):
    for name in DIST_COLUMN_PRIORITY:
        if name in df.columns:
            return name
    cand = [c for c in df.columns if c.lower().startswith("dist")]
    return cand[0] if cand else None

def compute_contact_weights(contacts_csv: Optional[str],
                            mode: str = "degree",
                            filter_min_seq_sep: int = 0,
                            max_dist: Optional[float] = None,
                            msa_length: Optional[int] = None):
    """
    Computes pair weights based on structural/evolutionary metrics.
    If mode='random', returns an empty set to signal the runner
    to generate null-control pairs.
    """
    if mode == "random":
        # Random mode doesn't use the CSV; it will be handled in the sampler
        return None, None

    if contacts_csv is None or not os.path.exists(contacts_csv):
        raise ValueError(f"CSV path required for mode '{mode}'")

    df = pd.read_csv(contacts_csv)

    # Standardize column names for residues
    cand_pairs = [("resi1", "resi2"), ("res_i", "res_j"), ("i", "j"), ("pos_i", "pos_j")]
    pair_cols = next(((a, b) for a, b in cand_pairs if a in df.columns and b in df.columns), None)

    if pair_cols is None:
        raise ValueError("Could not find residue index columns in CSV.")

    a_col, b_col = pair_cols
    df[a_col], df[b_col] = df[a_col].astype(int), df[b_col].astype(int)

    # Symmetry normalization (i < j)
    df = df.assign(_i=df[[a_col, b_col]].min(axis=1), _j=df[[a_col, b_col]].max(axis=1))

    # Apply Filters
    if filter_min_seq_sep > 0:
        df = df[(df["_j"] - df["_i"]).abs() >= filter_min_seq_sep]

    dcol = _pick_distance_column(df)
    if max_dist is not None:
        if dcol is None: raise ValueError("Distance filter requested but no 'dist' column found.")
        df = df[df[dcol].astype(float) <= float(max_dist)]

    df = df.drop_duplicates(subset=["_i", "_j"]).reset_index(drop=True)
    pairs = list(zip(df["_i"].tolist(), df["_j"].tolist()))

    if not pairs:
        raise ValueError("No contact pairs remain after filtering.")

    # Calculate Weights
    if mode == "uniform":
        w = np.ones(len(pairs))
    elif mode == "degree":
        deg = Counter()
        for p in pairs: deg.update(p)
        w = np.array([deg[i] + deg[j] for i, j in pairs], dtype=float)
    elif mode in ("dist", "hybrid"):
        invd = 1.0 / (df[dcol].astype(float).values + 1e-6)
        if mode == "dist":
            w = invd
        else: # hybrid
            deg = Counter()
            for p in pairs: deg.update(p)
            degw = np.array([deg[i] + deg[j] for i, j in pairs], dtype=float)
            w = degw * invd
    else:
        raise ValueError(f"Unknown mode: {mode}")

    return pairs, w / w.sum()

def weighted_sample_without_replacement(pairs, weights, k, rng_np):
    if k >= len(pairs):
        return pairs
    idxs = np.arange(len(pairs))
    chosen = rng_np.choice(idxs, size=k, replace=False, p=weights)
    return [pairs[i] for i in chosen]

# ---------- apply contacts ----------
def apply_contacts_in_msa(headers, seqs, contacts,
                          ref_index=0, one_based=True,
                          strength=0.6, seed=1234,
                          action="promote",
                          alpha=0.1, power_gamma=3.0,
                          overlap="sequential"):
    """
    overlap="sequential" (original behaviour): pairs are applied in order on shared
        column lists; a column shared by several pairs keeps the last pair's edit.
    overlap="lock" (review fix, opt-in): a column written by one pair is locked;
        later pairs that would write to a locked column are skipped (skipped=True).
    """
    if overlap not in ("sequential", "lock"):
        raise ValueError("overlap must be 'sequential' or 'lock'")
    locked = set()
    rng_py = random.Random(seed)
    rng_np = np.random.default_rng(seed)

    ref = seqs[ref_index]
    nogap2col = build_ref_nogap_to_col(ref, one_based)
    cols = [list(col) for col in zip(*seqs)]
    report = []

    for i_res, j_res in contacts:
        if i_res not in nogap2col or j_res not in nogap2col:
            continue
        ci, cj = nogap2col[i_res], nogap2col[j_res]
        colA, colB = cols[ci], cols[cj]

        # columnas que esta acción escribe (promote/break: B; promote_joint: A y B)
        written = {ci, cj} if action == "promote_joint" else {cj}
        if overlap == "lock" and (written & locked):
            mi = column_MI(colA, colB)
            report.append({"pair": (i_res, j_res), "mi_before": mi, "mi_after": mi,
                           "delta_mi": 0.0, "skipped": True})
            continue

        mi_before = column_MI(colA, colB)
        if action == "promote":
            promote_pair_columns(colA, colB, strength, rng_py, rng_np, alpha, power_gamma, ref_index)
        elif action == "break":
            break_pair_columns(colA, colB, strength, rng_py, ref_index)
        elif action == "promote_joint":
            #promote_by_joint(colA, colB, strength, rng_np, ref_index)
            promote_by_joint(colA, colB, strength, rng_np, ref_index, alpha, power_gamma)
        else:
            raise ValueError("action inválida.")

        mi_after = column_MI(colA, colB)
        if overlap == "lock":
            locked |= written
        report.append({
            "pair": (i_res, j_res),
            "mi_before": mi_before,
            "mi_after": mi_after,
            "delta_mi": mi_after - mi_before,
            "skipped": False
        })

    seqs_out = ["".join(row) for row in zip(*cols)]
    return headers, seqs_out, report

# ---------- experiment runner ----------
def make_unique_outdir(base_out_dir, msa_in_path, contacts_csv, action):
    msa_stem = os.path.splitext(os.path.basename(msa_in_path))[0]
    csv_stem = os.path.splitext(os.path.basename(contacts_csv))[0]
    out_dir = os.path.join(base_out_dir, f"{action}__{msa_stem}__{csv_stem}")
    os.makedirs(out_dir, exist_ok=True)
    return out_dir

def make_progressive_msas_weighted(msa_in_path, contacts_csv, base_out_dir,
                                   mode="hybrid",
                                   pair_fracs=(0.2, 0.4, 0.6, 0.8, 1.0),
                                   strength_grid=(0.6,),
                                   seeds=(11,),
                                   action="promote",
                                   ref_index=0, one_based=True,
                                   min_seq_sep: int = 6, # Default seq separation
                                   overlap: str = "sequential",
                                   write_reports: bool = True,
                                   **kwargs):

    headers, seqs = read_fasta_or_a3m(msa_in_path)
    headers, seqs = normalize_a3m(headers, seqs)

    # Calculate length of the reference sequence (excluding gaps)
    ref_seq = seqs[ref_index]
    nogap_indices = sorted(build_ref_nogap_to_col(ref_seq, one_based).keys())

    # 1. Get initial pairs/weights
    pairs, weights = compute_contact_weights(
        contacts_csv, mode=mode, filter_min_seq_sep=min_seq_sep,
        max_dist=kwargs.get("max_dist_for_gain") if kwargs.get("use_dist_filter") else None
    )

    # Determine reference count if we are in random mode
    if mode == "random":
        # To make it a fair control, we'll pick as many random pairs as exist in the original CSV
        ref_df = pd.read_csv(contacts_csv)
        total_pairs_to_sample_from = len(ref_df)
    else:
        total_pairs_to_sample_from = len(pairs)

    out_dir = make_unique_outdir(base_out_dir, msa_in_path, contacts_csv or "random", action)
    summary = []

    for frac in pair_fracs:
        k_pairs = max(1, int(round(frac * total_pairs_to_sample_from)))

        for strength in strength_grid:
            for seed in seeds:
                rng_np = np.random.default_rng(seed)
                rng_py = random.Random(seed)

                # --- Selection Logic ---
                if mode == "random":
                    # Generate totally random pairs that respect the min_seq_sep
                    chosen = []
                    while len(chosen) < k_pairs:
                        i, j = rng_np.choice(nogap_indices, size=2, replace=False)
                        if abs(i - j) >= min_seq_sep:
                            chosen.append((int(min(i, j)), int(max(i, j))))
                else:
                    chosen = weighted_sample_without_replacement(pairs, weights, k_pairs, rng_np)

                # Apply modifications
                h2, s2, rep = apply_contacts_in_msa(
                    headers, seqs, chosen, ref_index=ref_index,
                    one_based=one_based, strength=strength, seed=seed,
                    action=action, alpha=kwargs.get("alpha", 0.1),
                    power_gamma=kwargs.get("power_gamma", 3.0),
                    overlap=overlap
                )

                # --- Reporting & Saving ---
                tag = f"{mode}_frac{frac:.2f}_str{strength:.2f}_seed{seed}"
                write_fasta(os.path.join(out_dir, f"msa_{tag}.a3m"), h2, s2)
                if write_reports:
                    pd.DataFrame(rep).to_csv(os.path.join(out_dir, f"report_{tag}.csv"), index=False)

                dmi = [r["delta_mi"] for r in rep if not r.get("skipped")]
                summary.append(dict(
                    tag=tag, mode=mode, action=action, frac=frac,
                    strength=strength, seed=seed, k_pairs=len(chosen),
                    delta_mi_mean=float(np.mean(dmi)) if dmi else 0.0,
                    n_skipped=sum(1 for r in rep if r.get("skipped"))
                ))

    pd.DataFrame(summary).to_csv(os.path.join(out_dir, "summary.csv"), index=False)
    return out_dir


def _fmt(x):
    if isinstance(x, float):
        s = f"{x:.3g}"
    else:
        s = str(x)
    return s.replace(".", "_")

def copy_a3m_tree_to_master(src_root, dst_root, extra_suffix):
    os.makedirs(dst_root, exist_ok=True)
    for root, _, files in os.walk(src_root):
        for f in files:
            if not f.endswith(".a3m"):
                continue
            name, ext = os.path.splitext(f)
            dst = os.path.join(dst_root, f"{name}{extra_suffix}{ext}")
            shutil.copy2(os.path.join(root, f), dst)
