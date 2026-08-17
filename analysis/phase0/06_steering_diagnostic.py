"""
06_steering_diagnostic.py

Three tests to answer: "is the existing data strong enough to identify
column-level steering handles?"

  (1) Seed reproducibility:  for each (frac, str, A, G, MINSEP) cell,
      compute the variance in TM_Closed across the 3 seeds. If small,
      parameters dominate; if large, the modifier's randomness dominates.

  (2) Lasso TM_Closed ~ per-column edit fraction, inside the bimodal
      ensemble B_minus_A_control_negativo. CV-R^2 tells us how
      predictable the outcome is from MSA edits.

  (3) Extreme-pair contrast:  pick the most-closed and the most-open
      MSAs in B_minus_A_control_negativo and identify which columns
      have the biggest edit-fraction delta between them.
"""

import os, numpy as np, pandas as pd

# --- repo-relative paths -------------------------------------------------
from pathlib import Path as _Path
import os as _os
REPO   = _Path(__file__).resolve().parents[2]
RUNS   = _Path(_os.environ.get("ADK_RUNS",   REPO.parent / "nuevos_MSA" / "fuente_MSA_original" / "results"))
INPUTS = _Path(_os.environ.get("ADK_INPUTS", REPO / "data" / "AFF2024"))
# -------------------------------------------------------------------------


# --- mini ridge regression with k-fold CV (sklearn unavailable) ---
def ridge_fit(X, y, lam):
    n, p = X.shape
    A = X.T @ X + lam * np.eye(p)
    b = X.T @ y
    return np.linalg.solve(A, b)

def cv_r2(X, y, lam, k=5, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(y))
    folds = np.array_split(idx, k)
    r2s = []
    for i in range(k):
        te = folds[i]
        tr = np.concatenate([folds[j] for j in range(k) if j != i])
        Xtr, Xte = X[tr], X[te]
        ytr, yte = y[tr], y[te]
        mu = ytr.mean()
        beta = ridge_fit(Xtr, ytr - mu, lam)
        pred = Xte @ beta + mu
        ss_res = ((yte - pred)**2).sum()
        ss_tot = ((yte - yte.mean())**2).sum()
        r2s.append(1 - ss_res/ss_tot)
    return float(np.mean(r2s)), float(np.std(r2s))

def best_lambda(X, y, lam_grid=(0.01, 0.1, 1, 10, 100, 1000)):
    best = (-1e9, None)
    for lam in lam_grid:
        r2, _ = cv_r2(X, y, lam)
        if r2 > best[0]: best = (r2, lam)
    return best

def greedy_forward(X, y, max_features=15):
    """Forward stepwise selection by univariate R² of residual."""
    n, p = X.shape
    selected = []
    residual = y - y.mean()
    chosen_path = []
    for step in range(max_features):
        best_i, best_r = None, -1e9
        for j in range(p):
            if j in selected: continue
            xj = X[:, j]
            if xj.std() < 1e-9: continue
            # simple r²
            xj_c = xj - xj.mean()
            r = (xj_c @ residual) / (np.linalg.norm(xj_c) * np.linalg.norm(residual) + 1e-12)
            r2 = r**2
            if r2 > best_r: best_r, best_i = r2, j
        if best_i is None: break
        selected.append(best_i)
        # fit OLS on selected columns; recompute residual
        Xs = X[:, selected]
        beta = np.linalg.lstsq(Xs - Xs.mean(0), y - y.mean(), rcond=None)[0]
        residual = (y - y.mean()) - (Xs - Xs.mean(0)) @ beta
        cum_r2_train = 1 - (residual**2).sum() / ((y - y.mean())**2).sum()
        # CV r2 of the current selected set under ridge
        cv_r2_val, _ = cv_r2(Xs, y, lam=0.1)
        chosen_path.append((step+1, best_i, cum_r2_train, cv_r2_val))
    return selected, chosen_path

OUT  = str(REPO / "results" / "phase0")
M    = np.load(os.path.join(OUT,"col_edit_matrix.npy"))      # (1296, 214)
diff = pd.read_csv(os.path.join(OUT,"msa_diff_summary.csv"))
pred = pd.read_csv(os.path.join(OUT,"all_predictions_merged.csv"))

# domain helper from paper centroids
def domain_paper(res):
    if 35 <= res <= 55:   return "NMP"
    if 90 <= res <= 100:  return "CORE_anchor"
    if 115 <= res <= 124: return "hinge"
    if 125 <= res <= 153: return "LID"
    if 179 <= res <= 185: return "CORE2"
    return "CORE"

# join MSAs and predictions for B_minus_A_control specifically
pred["fname_base"] = pred["fname"].str.replace(".cif","",regex=False)
diff["fname_base"] = diff["fname"].str.replace(".a3m","",regex=False)
merged = diff.merge(pred[["fname_base","ensemble","NMP_Angle","LID_Angle",
                          "TM_Closed","TM_Open","dTM"]],
                    on=["fname_base","ensemble"], how="left")

# ==================================================================
# (1) SEED REPRODUCIBILITY
# ==================================================================
print("=" * 70)
print("(1) SEED REPRODUCIBILITY")
print("=" * 70)

# group by everything EXCEPT seed; we should have 3 rows per group
group_cols = ["ensemble","act","frac","strg","alpha","gamma","minsep"]
groups = merged.groupby(group_cols)
print(f"Number of (param,ensemble) cells: {groups.ngroups}")
print(f"Mean rows per cell: {groups.size().mean():.1f}  (should be 3 if balanced)")

cell_stats = groups.agg(
    n=("seed","count"),
    TM_C_mean=("TM_Closed","mean"),
    TM_C_std =("TM_Closed","std"),
    TM_O_std =("TM_Open","std"),
    NMP_std  =("NMP_Angle","std"),
    LID_std  =("LID_Angle","std"),
).reset_index()
# overall std (across all MSAs) for comparison
overall = merged.groupby("ensemble")[["TM_Closed","TM_Open","NMP_Angle","LID_Angle"]].std().add_prefix("overall_std_").reset_index()
cell_stats = cell_stats.merge(overall, on="ensemble")

# ratio of within-cell (seed-only) std to overall ensemble std
# small = seed barely matters; large = seed dominates
cell_stats["ratio_TM_C"]  = cell_stats["TM_C_std"] / cell_stats["overall_std_TM_Closed"]
cell_stats["ratio_NMP"]   = cell_stats["NMP_std"]  / cell_stats["overall_std_NMP_Angle"]
cell_stats["ratio_LID"]   = cell_stats["LID_std"]  / cell_stats["overall_std_LID_Angle"]

print("\nMedian seed-only / overall std ratio per ensemble:")
print(cell_stats.groupby("ensemble")[["ratio_TM_C","ratio_NMP","ratio_LID"]].median().round(3).to_string())
print("\n(ratio = 0 → seed doesn't matter; ratio = 1 → seed explains everything)")

# Inside B_minus_A_control specifically
bmac = cell_stats[cell_stats["ensemble"]=="B_minus_A_control_negativo"]
print("\nB_minus_A_control_negativo — distribution of within-seed TM_Closed std:")
print(bmac["TM_C_std"].describe().round(3).to_string())
print("Overall TM_Closed std in this ensemble:",
      merged[merged["ensemble"]=="B_minus_A_control_negativo"]["TM_Closed"].std().round(3))

# ==================================================================
# (2) LASSO inside B_minus_A_control_negativo
# ==================================================================
print("\n" + "=" * 70)
print("(2) LASSO  TM_Closed ~ per-column-edit-fraction  (B_minus_A_control)")
print("=" * 70)

sub = merged[merged["ensemble"]=="B_minus_A_control_negativo"]
idx = sub["mat_idx"].values
X = M[idx]                        # (324, 214)
y_closed = sub["TM_Closed"].values
y_open   = sub["TM_Open"].values
y_nmp    = sub["NMP_Angle"].values
y_lid    = sub["LID_Angle"].values

def predict_report(X, y, label, max_features=15):
    # standardise
    mu, sd = X.mean(0), X.std(0)
    sd_safe = np.where(sd > 1e-9, sd, 1.0)
    Xs = (X - mu) / sd_safe
    keep = sd > 1e-9
    Xs = Xs[:, keep]
    keep_idx = np.where(keep)[0]
    # CV R² with the full feature set (ridge)
    full_r2, full_sd = cv_r2(Xs, y, lam=10.0)
    # Best lambda
    best_r2, lam = best_lambda(Xs, y)
    # Greedy forward selection
    sel, path = greedy_forward(Xs, y, max_features=max_features)
    print(f"\n=== {label} ===")
    print(f"  ridge CV R² (all 214 columns)   = {full_r2:+.3f} ± {full_sd:.3f}")
    print(f"  ridge CV R² (best λ={lam})        = {best_r2:+.3f}")
    print(f"  greedy stepwise (resid → col → cum train R² → CV R²):")
    for step, j, train_r2, cv_r2_val in path:
        real_col = keep_idx[j]
        print(f"     step {step:2d}: res {real_col+1:3d} ({domain_paper(real_col+1):11s})  "
              f"cum train R²={train_r2:.3f}  CV R²={cv_r2_val:+.3f}")
    return path

predict_report(X, y_closed, "TM_Closed", max_features=15)
predict_report(X, y_open,   "TM_Open",   max_features=15)
predict_report(X, y_nmp,    "NMP_Angle", max_features=10)
predict_report(X, y_lid,    "LID_Angle", max_features=10)
predict_report(X, np.abs(sub["dTM"].values), "|dTM|", max_features=10)

# ==================================================================
# (3) EXTREME PAIR within B_minus_A_control
# ==================================================================
print("\n" + "=" * 70)
print("(3) EXTREME PAIR — most-closed vs most-open MSA inside B_minus_A_control")
print("=" * 70)

sub = sub.reset_index(drop=True)
i_closed = sub["TM_Closed"].idxmax()
i_open   = sub["TM_Open"].idxmax()
print(f"\nMost-closed MSA: {sub.loc[i_closed,'fname']}")
print(f"  TM_Closed={sub.loc[i_closed,'TM_Closed']:.3f}  TM_Open={sub.loc[i_closed,'TM_Open']:.3f}")
print(f"  params: frac={sub.loc[i_closed,'frac']}  str={sub.loc[i_closed,'strg']}  seed={sub.loc[i_closed,'seed']}  a={sub.loc[i_closed,'alpha']}  g={sub.loc[i_closed,'gamma']}  ms={sub.loc[i_closed,'minsep']}")
print(f"\nMost-open MSA: {sub.loc[i_open,'fname']}")
print(f"  TM_Closed={sub.loc[i_open,'TM_Closed']:.3f}  TM_Open={sub.loc[i_open,'TM_Open']:.3f}")
print(f"  params: frac={sub.loc[i_open,'frac']}  str={sub.loc[i_open,'strg']}  seed={sub.loc[i_open,'seed']}  a={sub.loc[i_open,'alpha']}  g={sub.loc[i_open,'gamma']}  ms={sub.loc[i_open,'minsep']}")

mat_c = sub.loc[i_closed,"mat_idx"]
mat_o = sub.loc[i_open,"mat_idx"]
ec, eo = M[mat_c], M[mat_o]
delta = ec - eo
print("\nTop 12 columns where 'closed' MSA edited MORE than 'open' MSA:")
for c in np.argsort(delta)[::-1][:12]:
    print(f"  res {c+1:3d} ({domain_paper(c+1):11s})  closed_ef={ec[c]:.3f}  open_ef={eo[c]:.3f}  Δ={delta[c]:+.3f}")
print("\nTop 12 columns where 'open' MSA edited MORE than 'closed' MSA:")
for c in np.argsort(delta)[:12]:
    print(f"  res {c+1:3d} ({domain_paper(c+1):11s})  closed_ef={ec[c]:.3f}  open_ef={eo[c]:.3f}  Δ={delta[c]:+.3f}")

# ==================================================================
# (3b) BROADER comparison: top-decile closed vs top-decile open
# ==================================================================
print("\n" + "=" * 70)
print("(3b) Mean edit-frac in top-10% closed vs top-10% open (B_minus_A_control)")
print("=" * 70)

q_closed = sub["TM_Closed"].quantile(0.9)
q_open   = sub["TM_Open"].quantile(0.9)
top_closed = sub[sub["TM_Closed"] >= q_closed]
top_open   = sub[sub["TM_Open"]   >= q_open]
print(f"n_top_closed={len(top_closed)}  n_top_open={len(top_open)}")

mc = M[top_closed["mat_idx"].values].mean(axis=0)
mo = M[top_open["mat_idx"].values].mean(axis=0)
delta = mc - mo
print("\nTop 12 columns enriched in the closed top-decile relative to open top-decile:")
for c in np.argsort(delta)[::-1][:12]:
    print(f"  res {c+1:3d} ({domain_paper(c+1):11s})  closed_top_ef={mc[c]:.3f}  open_top_ef={mo[c]:.3f}  Δ={delta[c]:+.3f}")
print("\nTop 12 columns enriched in the open top-decile relative to closed top-decile:")
for c in np.argsort(delta)[:12]:
    print(f"  res {c+1:3d} ({domain_paper(c+1):11s})  closed_top_ef={mc[c]:.3f}  open_top_ef={mo[c]:.3f}  Δ={delta[c]:+.3f}")
