"""
05_verify_and_export.py
- Verification pass that recomputes headline numbers independently.
- Re-tags every column with the paper's canonical AdK domain.
- Exports the final per-column tendency table that Maia can use to drive
  her MSA-modifier algorithm.
"""
import os, numpy as np, pandas as pd

# --- repo-relative paths -------------------------------------------------
from pathlib import Path as _Path
import os as _os
REPO   = _Path(__file__).resolve().parents[2]
RUNS   = _Path(_os.environ.get("ADK_RUNS",   REPO.parent / "nuevos_MSA" / "fuente_MSA_original" / "results"))
INPUTS = _Path(_os.environ.get("ADK_INPUTS", REPO / "data" / "AFF2024"))
# -------------------------------------------------------------------------


OUT  = str(REPO / "results" / "phase0")
M    = np.load(os.path.join(OUT,"col_edit_matrix.npy"))
diff = pd.read_csv(os.path.join(OUT,"msa_diff_summary.csv"))
pred = pd.read_csv(os.path.join(OUT,"all_predictions_merged.csv"))

# ---- domain map using the paper's centroids (1-based) ----
def domain_paper(res):
    if 35 <= res <= 55:   return "NMP"
    if 90 <= res <= 100:  return "CORE_anchor"
    if 115 <= res <= 124: return "hinge"          # corelid
    if 125 <= res <= 153: return "LID"
    if 179 <= res <= 185: return "CORE2"
    return "CORE"

# =================== verification ===================
print("== VERIFICATION ==")
# 1. recompute spread (NMP_std + LID_std) directly from pred
v_spread = pred.groupby("ensemble").apply(
    lambda d: d["NMP_Angle"].std() + d["LID_Angle"].std()
).rename("spread_idx_v2")
print("\nSpread index recomputed (should match 01_param_sweep):")
print(v_spread.to_string())

# 2. independently compute fraction closed using clearer thresholds
def state2(row):
    tc, to = row["TM_Closed"], row["TM_Open"]
    diff_  = tc - to
    if diff_ >  0.10:  return "closed"
    if diff_ < -0.10:  return "open"
    return "intermediate"
pred["state2"] = pred.apply(state2, axis=1)
print("\nState fractions (alt definition, |dTM|>0.10):")
print(pred.groupby("ensemble")["state2"].value_counts(normalize=True).unstack().round(3).to_string())

# 3. sanity check: query row (index 0) is never edited
print("\nMax query-row edits across all 1296 MSAs (should be 0):",
      diff["query_row_edits"].max())

# 4. sanity check: total edits asymmetry
print("\nMean total edits per ensemble:")
print(diff.groupby("ensemble")["total_edits"].mean().round(0).to_string())

# =================== export ===================
# Build a clean per-(ensemble, residue) tendency table
ens_list = pred["ensemble"].unique().tolist()
rows = []
for ens in ens_list:
    idx = diff[diff["ensemble"]==ens]["mat_idx"].values
    X   = M[idx]                                # (n, L)
    y_tc = pred.iloc[idx]["TM_Closed"].values
    y_to = pred.iloc[idx]["TM_Open"].values
    y_nmp= pred.iloc[idx]["NMP_Angle"].values
    y_lid= pred.iloc[idx]["LID_Angle"].values
    def pcorr(X, y):
        Xc = X - X.mean(0); yc = y - y.mean()
        den = np.sqrt((Xc**2).sum(0)) * np.sqrt((yc**2).sum())
        with np.errstate(invalid='ignore'):
            r = (Xc * yc[:,None]).sum(0) / den
        return np.nan_to_num(r, nan=0.0)
    r_tc = pcorr(X, y_tc); r_to = pcorr(X, y_to)
    r_nmp= pcorr(X, y_nmp); r_lid= pcorr(X, y_lid)
    for c in range(X.shape[1]):
        rows.append(dict(
            ensemble=ens, residue=c+1, domain=domain_paper(c+1),
            mean_edit_frac=X[:,c].mean(),
            r_TM_Closed=r_tc[c], r_TM_Open=r_to[c],
            r_NMP=r_nmp[c], r_LID=r_lid[c],
        ))
tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT,"final_tendency_table.csv"), index=False)
print("\nWrote:", os.path.join(OUT,"final_tendency_table.csv"))

# Concise "shortlist": for each ensemble, top edited columns AND top
# closed-promoter columns (high r_TM_Closed) AND top open-promoter columns.
print("\n=== SHORTLIST per ensemble ===")
for ens in ens_list:
    sub = tab[tab["ensemble"]==ens]
    print(f"\n--- {ens} ---")
    top_edit = sub.sort_values("mean_edit_frac", ascending=False).head(8)
    print("  most-edited columns:")
    for _, r in top_edit.iterrows():
        print(f"    res {int(r.residue):3d} ({r.domain:11s})  edit_frac={r.mean_edit_frac:.3f}")
    top_closed = sub.sort_values("r_TM_Closed", ascending=False).head(5)
    print("  edits at these cols correlate with CLOSED:")
    for _, r in top_closed.iterrows():
        print(f"    res {int(r.residue):3d} ({r.domain:11s})  r={r.r_TM_Closed:+.3f}")
    top_open = sub.sort_values("r_TM_Open", ascending=False).head(5)
    print("  edits at these cols correlate with OPEN:")
    for _, r in top_open.iterrows():
        print(f"    res {int(r.residue):3d} ({r.domain:11s})  r={r.r_TM_Open:+.3f}")

# per-domain reweighted summary
agg = tab.groupby(["ensemble","domain"])["mean_edit_frac"].mean().unstack().round(4)
print("\n=== Per-domain mean edit fraction (paper domains) ===")
print(agg.to_string())
agg.to_csv(os.path.join(OUT,"per_domain_summary_paper.csv"))
