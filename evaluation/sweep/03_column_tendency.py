"""
03_column_tendency.py

Given the per-MSA per-column edit-fraction matrix from 02_msa_diff.py
and the structural outcomes from 01_param_sweep.py, find:

  (a) Which COLUMNS are edited preferentially by each ensemble.
  (b) Which COLUMNS correlate (across MSAs) with high spread, with
      closed-state production, or with open-state production.
  (c) Where those columns sit in the protein domains.

AdK domain definition for chain A (residues 1-214) - based on Bhabha/Schultz/
Petsko literature & common AdK partitioning (1-based residue numbers):
  CORE    : 1-29, 60-121, 161-214
  NMP_bind: 30-59
  LID     : 122-159 (the mobile lid)
We'll annotate column index by adding 1 for residue number (query has no gaps).
"""

import os, numpy as np, pandas as pd

# --- repo-relative paths -------------------------------------------------
from pathlib import Path as _Path
import os as _os
REPO   = _Path(__file__).resolve().parents[2]
RUNS   = _Path(_os.environ.get("ADK_RUNS",   REPO.parent / "nuevos_MSA" / "fuente_MSA_original" / "results"))
INPUTS = _Path(_os.environ.get("ADK_INPUTS", REPO / "data" / "AFF2024"))
# -------------------------------------------------------------------------


OUT  = str(REPO / "results" / "sweep_evaluation")
M    = np.load(os.path.join(OUT, "col_edit_matrix.npy"))      # (1296, 214)
diff = pd.read_csv(os.path.join(OUT, "msa_diff_summary.csv"))
pred = pd.read_csv(os.path.join(OUT, "all_predictions_merged.csv"))

# join: msa_diff_summary.fname has .a3m; predictions has .cif → strip
pred["fname_base"] = pred["fname"].str.replace(".cif","",regex=False)
diff["fname_base"] = diff["fname"].str.replace(".a3m","",regex=False)
merged = diff.merge(pred[["fname_base","ensemble","NMP_Angle","LID_Angle","TM_Closed","TM_Open","dTM","state"]],
                    on=["fname_base","ensemble"], how="left")
print("Joined:", merged.shape, "missing matches:", merged["NMP_Angle"].isna().sum())

# domain map
def domain_of(res):
    if 1 <= res <= 29 or 60 <= res <= 121 or 161 <= res <= 214: return "CORE"
    if 30 <= res <= 59:   return "NMP"
    if 122 <= res <= 159: return "LID"
    return "?"
domains = np.array([domain_of(i+1) for i in range(M.shape[1])])

# ============================================================
# (a) Per-ensemble mean column-edit profile
# ============================================================
ens_profiles = {}
for ens in merged["ensemble"].unique():
    idx = merged[merged["ensemble"]==ens]["mat_idx"].values
    ens_profiles[ens] = M[idx].mean(axis=0)

prof_df = pd.DataFrame(ens_profiles)
prof_df["residue"] = np.arange(1, M.shape[1]+1)
prof_df["domain"]  = domains
prof_df = prof_df[["residue","domain"] + list(ens_profiles.keys())]
prof_df.to_csv(os.path.join(OUT,"per_column_profile.csv"), index=False)

print("\n--- Top-10 edited columns per ensemble (mean edit fraction) ---")
for ens, prof in ens_profiles.items():
    top = np.argsort(prof)[::-1][:10]
    print(f"\n[{ens}]")
    for c in top:
        print(f"   res {c+1:3d} ({domain_of(c+1)})  mean_edit_frac={prof[c]:.3f}")

# ============================================================
# (b) Per-column correlation with downstream outcomes
# ============================================================
# For each column, correlate edit_frac (across MSAs) with NMP_Angle,
# LID_Angle, TM_Closed, TM_Open, |dTM|, state-counts.
#
# We do this WITHIN each ensemble (so the "act" doesn't bleed in).
# ============================================================

corr_rows = []
for ens in merged["ensemble"].unique():
    sub = merged[merged["ensemble"]==ens]
    idx = sub["mat_idx"].values
    X = M[idx]                              # (n_msa, L)
    for outcome in ["NMP_Angle","LID_Angle","TM_Closed","TM_Open","dTM"]:
        y = sub[outcome].values
        # Pearson per column
        Xc = X - X.mean(axis=0)
        yc = y - y.mean()
        denom = (np.sqrt((Xc**2).sum(axis=0)) * np.sqrt((yc**2).sum()))
        with np.errstate(invalid='ignore'):
            r = (Xc * yc[:,None]).sum(axis=0) / denom
        r = np.nan_to_num(r, nan=0.0)
        for c in range(M.shape[1]):
            corr_rows.append(dict(
                ensemble=ens, outcome=outcome,
                residue=c+1, domain=domain_of(c+1), r=r[c]
            ))
corr_df = pd.DataFrame(corr_rows)
corr_df.to_csv(os.path.join(OUT,"per_column_correlations.csv"), index=False)

# ============================================================
# Concise summary: per ensemble, top |r| columns for each outcome
# ============================================================
print("\n=== Per ensemble: columns whose edit-fraction MOST correlates with each outcome ===")
for ens in corr_df["ensemble"].unique():
    print(f"\n--- {ens} ---")
    for outc in ["NMP_Angle","LID_Angle","TM_Closed","TM_Open","dTM"]:
        s = corr_df[(corr_df["ensemble"]==ens) & (corr_df["outcome"]==outc)]
        s = s.reindex(s["r"].abs().sort_values(ascending=False).index).head(5)
        s = s[s["r"].abs() > 0.05]
        if s.empty: continue
        print(f"  [{outc}]")
        for _, row in s.iterrows():
            print(f"    res {int(row['residue']):3d} ({row['domain']})  r={row['r']:+.3f}")

# ============================================================
# Per-domain aggregate: total mean edits per ensemble per domain
# ============================================================
print("\n=== Per-domain mean column-edit fraction by ensemble ===")
dom_summary = []
for ens, prof in ens_profiles.items():
    row = {"ensemble":ens}
    for d in ["NMP","LID","CORE"]:
        mask = domains==d
        row[d+"_mean_edit_frac"] = prof[mask].mean()
        row[d+"_sum_edit_frac"]  = prof[mask].sum()
    dom_summary.append(row)
dom_df = pd.DataFrame(dom_summary)
print(dom_df.to_string(index=False))
dom_df.to_csv(os.path.join(OUT,"per_domain_summary.csv"), index=False)
