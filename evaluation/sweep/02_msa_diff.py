"""
02_msa_diff.py
Diff every modified MSA against the original, producing a per-column
edit-fraction matrix (1296 MSAs × 214 columns) plus per-MSA summary
stats. This is the substrate for the tendency analysis.
"""

import os, re, glob, json
import numpy as np
import pandas as pd

# --- repo-relative paths -------------------------------------------------
from pathlib import Path as _Path
import os as _os
REPO   = _Path(__file__).resolve().parents[2]
RUNS   = _Path(_os.environ.get("ADK_RUNS",   REPO.parent / "nuevos_MSA" / "fuente_MSA_original" / "results"))
INPUTS = _Path(_os.environ.get("ADK_INPUTS", REPO / "data" / "AFF2024"))
# -------------------------------------------------------------------------


ROOT  = str(RUNS.parent)
ORIG  = str(INPUTS / "AdK" / "Frustration-Energy-space" / "MSA.a3m")
OUT   = str(REPO / "results" / "sweep_evaluation")

ENS = [
    "all_runs_a3m_A_minus_B",
    "all_runs_a3m_A_minus_B_control_negativo",
    "all_runs_a3m_B_minus_A",
    "all_runs_a3m_B_minus_A_control_negativo",
]

FNAME_RE = re.compile(
    r"msa_(?P<mode>\w+?)_frac(?P<frac>[0-9.]+)_str(?P<strg>[0-9.]+)_seed(?P<seed>\d+)"
    r"__ACT-(?P<act>[A-Za-z_]+)__MODE-(?P<mode2>\w+?)"
    r"__A-(?P<alpha>[0-9_]+?)__G-(?P<gamma>\d+)__MINSEP-(?P<minsep>\d+)\.a3m"
)

def load_msa_chars(path, L=214):
    """Fast loader. Returns np.array of shape (N, L) dtype '|S1'."""
    seqs = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i % 2 == 1:                       # sequence line
                seqs.append(line.rstrip())
    arr = np.frombuffer("".join(seqs).encode(), dtype="|S1").reshape(len(seqs), -1)
    return None, arr

def parse_fname(fn):
    m = FNAME_RE.match(os.path.basename(fn))
    if not m: return None
    g = m.groupdict()
    return dict(
        frac   = float(g["frac"]),
        strg   = float(g["strg"]),
        seed   = int(g["seed"]),
        act    = g["act"],
        alpha  = float(g["alpha"].replace("_", ".").rstrip(".")),
        gamma  = int(g["gamma"]),
        minsep = int(g["minsep"]),
    )

# ---- load original once ----
orig_headers, orig = load_msa_chars(ORIG)
N, L = orig.shape
print(f"Original MSA: {N} sequences × {L} columns")

# ---- iterate every modified MSA ----
col_edit_matrix = []          # (M, L) per-column edit FRACTION
row_edit_summary = []         # per-MSA aggregate
row_max_changed  = []         # which row got hit the hardest
sub_records       = []        # for top-substitution histogram (sampled)

records = []
for ens in ENS:
    files = sorted(glob.glob(os.path.join(ROOT, ens, "*.a3m")))
    print(f"  {ens}: {len(files)} files")
    for fp in files:
        _, mod = load_msa_chars(fp)
        if mod.shape != orig.shape:
            print(" ! shape mismatch", fp, mod.shape); continue
        diff = (mod != orig)                 # bool (N, L)
        col_edits = diff.sum(axis=0)         # (L,)  total per column
        row_edits = diff.sum(axis=1)         # (N,)  total per row
        total     = diff.sum()
        # gap introductions / removals
        gap_intro = ((mod == b"-") & (orig != b"-")).sum()
        gap_remov = ((mod != b"-") & (orig == b"-")).sum()
        params = parse_fname(fp) or {}
        rec = dict(
            ensemble = ens.replace("all_runs_a3m_",""),
            fname    = os.path.basename(fp),
            **params,
            total_edits      = int(total),
            edited_cells_pct = total/(N*L),
            rows_edited      = int((row_edits>0).sum()),
            cols_edited      = int((col_edits>0).sum()),
            max_col_edits    = int(col_edits.max()),
            argmax_col       = int(col_edits.argmax()),
            max_row_edits    = int(row_edits.max()),
            gap_intro        = int(gap_intro),
            gap_remov        = int(gap_remov),
            # query row (index 0) edits  ← KEY: did they touch the query?
            query_row_edits  = int(diff[0].sum()),
        )
        records.append(rec)
        col_edit_matrix.append(col_edits / N)   # fraction of rows edited at col

col_edit_matrix = np.array(col_edit_matrix)  # (M, L)
np.save(os.path.join(OUT, "col_edit_matrix.npy"), col_edit_matrix.astype("float32"))

df = pd.DataFrame(records)
df.to_csv(os.path.join(OUT, "msa_diff_summary.csv"), index=False)

# attach an ordinal index aligning the matrix rows with df rows
df["mat_idx"] = np.arange(len(df))
df.to_csv(os.path.join(OUT, "msa_diff_summary.csv"), index=False)

print("\nSaved:")
print(" ", os.path.join(OUT, "col_edit_matrix.npy"), col_edit_matrix.shape)
print(" ", os.path.join(OUT, "msa_diff_summary.csv"))
print("\nFirst rows:")
print(df.head().to_string(index=False))
print("\nAggregate per ensemble:")
print(df.groupby("ensemble")[["total_edits","cols_edited","rows_edited",
                              "max_col_edits","max_row_edits","query_row_edits",
                              "gap_intro","gap_remov"]].mean().round(2).to_string())
