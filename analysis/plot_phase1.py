#!/usr/bin/env python3
"""
plot_phase1.py — figures for the v2_phase1 hypothesis analysis.

Produces:
  fig_dose_response_dTM.png   — 4×7 panel: dose vs dTM, per (col × theory × dir)
  fig_dose_response_LID.png   — same, LID angle
  fig_mechanism_dose100.png   — A vs B vs C at dose 100, dTM and LID
  fig_seed_dominance.png      — distribution of dTM per Boltz seed across all 152 edits
  fig_summary_2D.png          — NMP vs LID, baseline + Phase-1, colored by theory
"""
from __future__ import annotations
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ANA  = str(REPO / "results" / "phase1")

df    = pd.read_csv(f"{ANA}/per_structure_metrics.csv")
cells = pd.read_csv(f"{ANA}/per_cell_stats.csv")
with open(f"{ANA}/baseline_reference.json") as fh:
    base = json.load(fh)
B = base["metrics"]

# Reference points (1AKE / 4AKE) — computed by score_boltz on the same refs
NMP_C, LID_C = 45.14, 106.40
NMP_O, LID_O = 73.77, 146.22

COLS = [50, 64, 136, 185, 190, 202, 206]
THEORIES = ["A", "B", "C"]
DIRS = ["closed", "open"]
DOSES = [25, 50, 75, 100]
THEORY_COLOR = {"A": "#1f77b4", "B": "#d62728", "C": "#2ca02c"}
DIR_LS = {"closed": "-", "open": "--"}


def baseline_band(ax, metric):
    m = B[metric]
    ax.axhline(m["mean"], color="grey", lw=0.8)
    ax.axhspan(m["mean"] - m["sd"], m["mean"] + m["sd"], alpha=0.10, color="grey")


# ---------------------------------------------------------------------------
# Figure 1: dose-response dTM
# ---------------------------------------------------------------------------
def fig_dose_response(metric: str, ylim, fname: str, title: str):
    fig, axs = plt.subplots(2, 7, figsize=(20, 7), sharey=True, sharex=True)
    for i, dr in enumerate(DIRS):
        for j, col in enumerate(COLS):
            ax = axs[i, j]
            baseline_band(ax, metric)
            for th in THEORIES:
                sub = cells[(cells.col == col) & (cells.theory == th) &
                            (cells.direction == dr)].sort_values("dose")
                if sub.empty:
                    continue
                y = sub[f"{metric}_mean"].values
                e = sub[f"{metric}_sd"].values / np.sqrt(5)  # SE for n=5
                ax.errorbar(sub["dose"], y, yerr=e,
                            marker="o", color=THEORY_COLOR[th], lw=1.8,
                            label=f"Theory {th}", capsize=2)
            ax.set_title(f"col {col}  dir={dr}", fontsize=9)
            ax.set_xticks(DOSES)
            ax.set_xlim(15, 110)
            if ylim:
                ax.set_ylim(*ylim)
            if j == 0:
                ax.set_ylabel(f"{metric} mean ± SE")
            if i == 1:
                ax.set_xlabel("dose (% editable rows)")
    axs[0, -1].legend(loc="upper right", fontsize=8, framealpha=0.9)
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{ANA}/{fname}", dpi=130)
    plt.close(fig)
    print(f"wrote {fname}")


fig_dose_response("dTM", (-0.10, 0.20),
                  "fig_dose_response_dTM.png",
                  "Dose-response per (column × theory × direction): dTM = TM_Open − TM_Closed "
                  "(grey = unedited-MSA baseline mean ± 1 SD)")

fig_dose_response("LID_Angle", (95, 150),
                  "fig_dose_response_LID.png",
                  "Dose-response per (column × theory × direction): LID-Core angle "
                  "(grey = baseline; 4AKE=146°, 1AKE=106°)")

# ---------------------------------------------------------------------------
# Figure 2: mechanism test — A vs B vs C at dose 100
# ---------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
for ax, met, ylab in zip(axs, ["dTM", "LID_Angle"], ["dTM (Open − Closed)", "LID-Core angle"]):
    baseline_band(ax, met)
    cells_d100 = cells[cells.dose == 100]
    width = 0.12
    x = np.arange(len(COLS))
    for i, th in enumerate(THEORIES):
        for k, dr in enumerate(DIRS):
            sub = cells_d100[(cells_d100.theory == th) &
                             (cells_d100.direction == dr)].set_index("col")
            offset = (i - 1) * 0.30 + (k - 0.5) * width
            vals = [sub.loc[c, f"{met}_mean"] if c in sub.index else np.nan for c in COLS]
            errs = [sub.loc[c, f"{met}_sd"] / np.sqrt(5) if c in sub.index else 0 for c in COLS]
            ax.errorbar(x + offset, vals, yerr=errs,
                        marker="o" if dr == "closed" else "s",
                        ls="None", color=THEORY_COLOR[th], capsize=2,
                        label=f"{th} / {dr}", markersize=6)
    ax.set_xticks(x)
    ax.set_xticklabels([f"col {c}" for c in COLS])
    ax.set_ylabel(ylab)
    if met == "LID_Angle":
        ax.axhline(LID_C, color="red", lw=0.8, alpha=0.6); ax.text(0.5, LID_C - 1, "1AKE (closed)", color="red", fontsize=8)
        ax.axhline(LID_O, color="red", lw=0.8, ls="--", alpha=0.6); ax.text(0.5, LID_O - 4, "4AKE (open)", color="red", fontsize=8)
    ax.set_title(f"Mechanism test at dose 100: {met}")
    ax.legend(fontsize=7, ncol=3, loc="upper right")
fig.tight_layout()
fig.savefig(f"{ANA}/fig_mechanism_dose100.png", dpi=130)
plt.close(fig)
print("wrote fig_mechanism_dose100.png")

# ---------------------------------------------------------------------------
# Figure 3: seed dominance
# ---------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(13, 4.5))
for ax, met, ylab in zip(axs, ["dTM", "LID_Angle"], ["dTM", "LID-Core angle"]):
    seeds = sorted(df.boltz_seed.unique())
    for s in seeds:
        vals = df[df.boltz_seed == s][met].values
        ax.violinplot([vals], positions=[s], showmeans=True, widths=0.7)
    baseline_band(ax, met)
    ax.set_xticks(seeds)
    ax.set_xlabel("Boltz seed")
    ax.set_ylabel(ylab)
    ax.set_title(f"{met} per Boltz seed, across all 152 modified MSAs\n"
                 f"(if edits mattered, distributions wouldn't be tight)")
fig.tight_layout()
fig.savefig(f"{ANA}/fig_seed_dominance.png", dpi=130)
plt.close(fig)
print("wrote fig_seed_dominance.png")

# ---------------------------------------------------------------------------
# Figure 4: NMP vs LID 2D scatter, baseline-vs-Phase1
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 7))
# Baseline: load the raw CSV
base_csv = str(REPO / "results" / "phase0" / "baseline_control_314.csv")
try:
    bdf = pd.read_csv(base_csv)
    ax.scatter(bdf["NMP_Angle"], bdf["LID_Angle"], s=12, alpha=0.35,
               color="grey", label=f"Unedited-MSA baseline (n={len(bdf)})")
except Exception as e:
    print("baseline CSV not available:", e)

for th, color in THEORY_COLOR.items():
    sub = df[df.theory == th]
    ax.scatter(sub["NMP_Angle"], sub["LID_Angle"], s=12, alpha=0.55,
               color=color, label=f"Theory {th} (n={len(sub)})")
# Refs
ax.scatter([NMP_C], [LID_C], s=220, c="red", marker="o",
           edgecolors="darkred", zorder=5, label="1AKE (closed)")
ax.scatter([NMP_O], [LID_O], s=220, facecolors="none", edgecolors="red",
           linewidths=2.0, zorder=5, label="4AKE (open)")
ax.set_xlabel("NMP-Core angle (°)")
ax.set_ylabel("LID-Core angle (°)")
ax.set_title("AdK conformational space — unedited-MSA baseline vs Phase-1 edits (all Boltz-2)")
ax.legend(loc="lower right", fontsize=9)
fig.tight_layout()
fig.savefig(f"{ANA}/fig_summary_2D.png", dpi=130)
plt.close(fig)
print("wrote fig_summary_2D.png")
