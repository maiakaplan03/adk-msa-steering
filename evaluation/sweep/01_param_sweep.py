"""
01_param_sweep.py
Diagnostic pass over the 4 ensemble result CSVs.

For each ensemble (A-B, A-B_control, B-A, B-A_control) we compute:
  - global SPREAD metrics  (stdev in NMP/LID/TM_Closed/TM_Open + 2D bbox area
    in NMP-LID space)
  - DIRECTIONALITY        (mean TM_Open - TM_Closed, plus class counts:
    closed-like / open-like / intermediate / spread)
  - per-parameter slice    (mean & std of the same metrics grouped by
    frac, str, A, G, MINSEP, seed)

Outputs are written next to this script as CSV so the next step can pick
the most informative subset of MSAs for the content tendency hunt.
"""

import os
import re
import numpy as np
import pandas as pd

# --- repo-relative paths -------------------------------------------------
from pathlib import Path as _Path
import os as _os
REPO   = _Path(__file__).resolve().parents[2]
RUNS   = _Path(_os.environ.get("ADK_RUNS",   REPO.parent / "nuevos_MSA" / "fuente_MSA_original" / "results"))
INPUTS = _Path(_os.environ.get("ADK_INPUTS", REPO / "data" / "AFF2024"))
# -------------------------------------------------------------------------


# --- paths ---
RES_ROOT = str(RUNS)
OUT_DIR  = str(REPO / "results" / "sweep_evaluation")
os.makedirs(OUT_DIR, exist_ok=True)

ENSEMBLES = [
    "all_runs_a3m_A_minus_B_results",
    "all_runs_a3m_A_minus_B_control_negativo_results",
    "all_runs_a3m_B_minus_A_results",
    "all_runs_a3m_B_minus_A_control_negativo_results",
]

# parse filename to recover sweep parameters that the CSV doesn't always carry
# example: msa_degree_frac0.30_str0.30_seed11__ACT-break__MODE-degree__A-0_05__G-2__MINSEP-12.cif
FNAME_RE = re.compile(
    r"msa_(?P<mode>\w+?)_frac(?P<frac>[0-9.]+)_str(?P<strg>[0-9.]+)_seed(?P<seed>\d+)"
    r"__ACT-(?P<act>[A-Za-z_]+)__MODE-(?P<mode2>\w+?)"
    r"__A-(?P<alpha>[0-9_]+?)__G-(?P<gamma>\d+)__MINSEP-(?P<minsep>\d+)\.cif"
)

def parse_fname(fn):
    m = FNAME_RE.match(fn)
    if not m:
        return None
    g = m.groupdict()
    return dict(
        mode=g["mode"],
        frac=float(g["frac"]),
        strg=float(g["strg"]),
        seed=int(g["seed"]),
        act=g["act"],
        alpha=float(g["alpha"].replace("_", ".").rstrip(".")),
        gamma=int(g["gamma"]),
        minsep=int(g["minsep"]),
    )

# --- load + tidy ---
frames = []
for ens in ENSEMBLES:
    p = os.path.join(RES_ROOT, ens, "interactive_data_backup.csv")
    df = pd.read_csv(p)
    df["ensemble"] = ens.replace("all_runs_a3m_", "").replace("_results", "")
    # re-parse so we have clean dtypes
    parsed = df["fname"].apply(parse_fname).apply(pd.Series)
    df = pd.concat([df.drop(columns=[c for c in ["frac","str","A","G","MINSEP"] if c in df.columns]),
                    parsed], axis=1)
    frames.append(df)

all_df = pd.concat(frames, ignore_index=True)
all_df["dTM"] = all_df["TM_Open"] - all_df["TM_Closed"]   # >0 → open-like

print("Total rows:", len(all_df))
print("Per ensemble:", all_df["ensemble"].value_counts().to_dict())
print("Per action:",   all_df["act"].value_counts().to_dict())

# --- structure classification (used for "spread" intuition) ---
# closed-like: TM_Closed >> TM_Open
# open-like:   TM_Open   >> TM_Closed
# intermediate: both moderate
# We pick simple thresholds; tweak later.
def classify(row, hi=0.85, lo=0.75):
    tc, to = row["TM_Closed"], row["TM_Open"]
    if tc >= hi and to <  lo:        return "closed"
    if to >= hi and tc <  lo:        return "open"
    if tc >= hi and to >= hi:        return "both"     # near both refs (ambiguous)
    return "intermediate"

all_df["state"] = all_df.apply(classify, axis=1)

# --- global summary per ensemble ---
def ensemble_summary(df):
    out = {}
    out["n"]            = len(df)
    out["NMP_mean"]     = df["NMP_Angle"].mean()
    out["NMP_std"]      = df["NMP_Angle"].std()
    out["LID_mean"]     = df["LID_Angle"].mean()
    out["LID_std"]      = df["LID_Angle"].std()
    out["TM_C_mean"]    = df["TM_Closed"].mean()
    out["TM_O_mean"]    = df["TM_Open"].mean()
    out["dTM_mean"]     = df["dTM"].mean()           # directionality
    out["dTM_std"]      = df["dTM"].std()
    # 2D spread in NMP-LID space
    bbox = (df["NMP_Angle"].max()-df["NMP_Angle"].min()) * \
           (df["LID_Angle"].max()-df["LID_Angle"].min())
    out["NMP_LID_bbox"] = bbox
    # convex-hull-ish "spread"  = sum of stds (cheap)
    out["spread_idx"]   = df["NMP_Angle"].std() + df["LID_Angle"].std()
    # state distribution
    for s in ["closed", "open", "intermediate", "both"]:
        out[f"frac_{s}"] = (df["state"] == s).mean()
    return pd.Series(out)

ens_summary = all_df.groupby("ensemble").apply(ensemble_summary).reset_index()
ens_summary.to_csv(os.path.join(OUT_DIR, "ensemble_summary.csv"), index=False)
print("\n=== Ensemble summary ===")
print(ens_summary.to_string(index=False))

# --- per-parameter slice (within each ensemble) ---
PARAMS = ["frac", "strg", "alpha", "gamma", "minsep", "seed"]
rows = []
for ens, sub in all_df.groupby("ensemble"):
    for p in PARAMS:
        for val, ssub in sub.groupby(p):
            rows.append(dict(
                ensemble = ens,
                param    = p,
                value    = val,
                n        = len(ssub),
                NMP_std  = ssub["NMP_Angle"].std(),
                LID_std  = ssub["LID_Angle"].std(),
                spread_idx = ssub["NMP_Angle"].std() + ssub["LID_Angle"].std(),
                dTM_mean = ssub["dTM"].mean(),
                dTM_std  = ssub["dTM"].std(),
                frac_closed = (ssub["state"]=="closed").mean(),
                frac_open   = (ssub["state"]=="open").mean(),
                frac_inter  = (ssub["state"]=="intermediate").mean(),
            ))
param_slice = pd.DataFrame(rows)
param_slice.to_csv(os.path.join(OUT_DIR, "param_slice.csv"), index=False)
print("\n=== Per-parameter slice (head) ===")
print(param_slice.head(20).to_string(index=False))

# --- save the cleaned, merged table for the next stages ---
all_df.to_csv(os.path.join(OUT_DIR, "all_predictions_merged.csv"), index=False)
print("\nWrote merged table:", os.path.join(OUT_DIR, "all_predictions_merged.csv"))
