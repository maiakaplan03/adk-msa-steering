"""
04_plots.py
Produce the figures we'll hand back to Maia.
"""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

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
prof = pd.read_csv(os.path.join(OUT,"per_column_profile.csv"))
corr = pd.read_csv(os.path.join(OUT,"per_column_correlations.csv"))

ORDER = ["A_minus_B","A_minus_B_control_negativo","B_minus_A","B_minus_A_control_negativo"]
COLORS = {"A_minus_B":"#1f77b4","A_minus_B_control_negativo":"#aec7e8",
          "B_minus_A":"#d62728","B_minus_A_control_negativo":"#ff9896"}

# Domain centroids from angle_analysis_plot_fixed.py (1-based residues here)
#   nmp     35–55     (NMP_bind)
#   core    90–100    (CORE, vertex of NMP angle)
#   corelid 115–125   (CORE/LID hinge)
#   lid     125–153   (LID body)
#   core2   179–185   (CORE C-terminus arm)
DOMAINS = [
    ("CORE",   1,  34, "#cfd8dc"),
    ("NMP",   35,  55, "#ffe082"),
    ("CORE",  56, 114, "#cfd8dc"),
    ("hinge",115, 124, "#fbc02d"),   # corelid hinge region
    ("LID",  125, 153, "#a5d6a7"),
    ("CORE",154, 178, "#cfd8dc"),
    ("CORE2",179, 185, "#90a4ae"),
    ("CORE",186, 214, "#cfd8dc"),
]

mpl.rcParams.update({"font.size":9, "axes.spines.top":False,"axes.spines.right":False})

# ============================================================
# Fig 1.  NMP vs LID hinge scatter per ensemble
# ============================================================
fig, ax = plt.subplots(1, 4, figsize=(15, 3.5), sharex=True, sharey=True)
for i, ens in enumerate(ORDER):
    sub = pred[pred["ensemble"]==ens]
    ax[i].scatter(sub["NMP_Angle"], sub["LID_Angle"], s=8, alpha=0.6, c=COLORS[ens])
    ax[i].set_title(f"{ens}\nn={len(sub)}  NMP σ={sub['NMP_Angle'].std():.1f}°  LID σ={sub['LID_Angle'].std():.1f}°", fontsize=9)
    ax[i].set_xlabel("NMP hinge angle (°)")
    if i==0: ax[i].set_ylabel("LID hinge angle (°)")
plt.tight_layout()
plt.savefig(os.path.join(OUT,"fig1_NMP_vs_LID.png"), dpi=150)
plt.close()

# ============================================================
# Fig 2.  Per-column edit profile (mean across MSAs) per ensemble
# ============================================================
fig, ax = plt.subplots(4, 1, figsize=(13, 7), sharex=True)
for i, ens in enumerate(ORDER):
    p = prof[ens].values
    ax[i].bar(prof["residue"], p, width=1.0, color=COLORS[ens], edgecolor="none")
    # domain shading
    for dname, s, e, col in DOMAINS:
        ax[i].axvspan(s-0.5, e+0.5, color=col, alpha=0.25, zorder=0)
    ax[i].set_ylabel(ens.replace("_minus_","−").replace("_negativo",""), fontsize=8)
    ax[i].set_ylim(0, prof[ORDER].values.max()*1.05)
    # annotate top 5
    top = np.argsort(p)[::-1][:5]
    for t in top:
        ax[i].annotate(f"{t+1}", (t+1, p[t]), fontsize=7, ha="center", va="bottom")
ax[-1].set_xlabel("Residue number (query column)")
ax[0].set_title("Mean per-column edit fraction across 324 MSAs / ensemble  (yellow=NMP, green=LID, grey=CORE)")
plt.tight_layout()
plt.savefig(os.path.join(OUT,"fig2_column_edit_profile.png"), dpi=150)
plt.close()

# ============================================================
# Fig 3.  Per-column |correlation| with TM_Closed / TM_Open (B_minus_A_control)
# ============================================================
fig, ax = plt.subplots(2, 1, figsize=(13, 5), sharex=True)
for j, outc in enumerate(["TM_Closed","TM_Open"]):
    s = corr[(corr["ensemble"]=="B_minus_A_control_negativo") & (corr["outcome"]==outc)]
    s = s.sort_values("residue")
    bar_col = np.where(s["r"]>0, "#2e7d32", "#c62828")  # green=positive (more edits→more of that state)
    ax[j].bar(s["residue"], s["r"], width=1.0, color=bar_col, edgecolor="none")
    for dname, st, e, col in DOMAINS:
        ax[j].axvspan(st-0.5, e+0.5, color=col, alpha=0.25, zorder=0)
    ax[j].axhline(0, color="k", lw=0.4)
    ax[j].set_ylabel(f"r(edit-frac, {outc})")
    # annotate top 5 |r|
    top5 = s.reindex(s["r"].abs().sort_values(ascending=False).index).head(5)
    for _, row in top5.iterrows():
        ax[j].annotate(f"{int(row['residue'])}", (row["residue"], row["r"]),
                       fontsize=7, ha="center",
                       va="bottom" if row["r"]>0 else "top")
ax[-1].set_xlabel("Residue number")
ax[0].set_title("B_minus_A_control_negativo: per-column Pearson r between edit fraction and TM score\n(green = edits at this column tend to PRODUCE that state)")
plt.tight_layout()
plt.savefig(os.path.join(OUT,"fig3_correlation_with_TM.png"), dpi=150)
plt.close()

# ============================================================
# Fig 4.  Per-domain mean edits — bar plot per ensemble
# ============================================================
dom_df = pd.read_csv(os.path.join(OUT,"per_domain_summary.csv"))
fig, ax = plt.subplots(figsize=(8,4))
x = np.arange(len(ORDER))
w = 0.27
for k, dom in enumerate(["NMP","LID","CORE"]):
    y = [dom_df[dom_df["ensemble"]==e][f"{dom}_mean_edit_frac"].values[0] for e in ORDER]
    ax.bar(x + (k-1)*w, y, width=w, label=dom,
           color={"NMP":"#ffe082","LID":"#a5d6a7","CORE":"#cfd8dc"}[dom],
           edgecolor="k")
ax.set_xticks(x)
ax.set_xticklabels([e.replace("_minus_","−").replace("_negativo","\n(break)") for e in ORDER], fontsize=8)
ax.set_ylabel("Mean edit fraction per column")
ax.legend(loc="upper right")
ax.set_title("Per-domain MSA-editing intensity by ensemble\n(B−A_control suppresses NMP edits → conformational freedom)")
plt.tight_layout()
plt.savefig(os.path.join(OUT,"fig4_domain_intensity.png"), dpi=150)
plt.close()

print("Saved:")
for f in ["fig1_NMP_vs_LID.png","fig2_column_edit_profile.png",
          "fig3_correlation_with_TM.png","fig4_domain_intensity.png"]:
    print(" ", os.path.join(OUT,f))
