#!/usr/bin/env python3
"""
mi_on_real_msa.py — does promote_joint raise MI on the real AdK alignment, as it does on
the toy alignment of pipeline/03_sanity_check_mi.py?

Reads the per-run summary.csv files that pipeline/04_generate_sweep.py keeps with
--keep-runs, and tabulates the mean ΔMI over edited pairs for every configuration.

  python pipeline/04_generate_sweep.py --msa <MSA.a3m> --ensemble A_minus_B --out-root /tmp/mi --keep-runs
  python pipeline/04_generate_sweep.py --msa <MSA.a3m> --ensemble B_minus_A --out-root /tmp/mi --keep-runs
  python evaluation/mi_on_real_msa.py --runs-root /tmp/mi

Writes results/sanity_check/mi_on_real_msa.csv (one row per alignment).
"""
import argparse, glob, os, re
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument("--runs-root", required=True)
ap.add_argument("--out", default=os.path.join(REPO, "results", "sanity_check", "mi_on_real_msa.csv"))
args = ap.parse_args()

pat = re.compile(r"all_runs_a3m_(\w+?)__runs/run__act-(\w+?)__mode-(\w+?)__a-([\d_]+)__g-([\d_]+)__minsep-(\d+)")
rows = []
for f in sorted(glob.glob(os.path.join(args.runs_root, "*__runs", "*", "*", "summary.csv"))):
    m = pat.search(f)
    d = pd.read_csv(f)
    d["ensemble"] = m.group(1)
    d["alpha"] = float(m.group(4).replace("_", "."))
    d["gamma"] = float(m.group(5).replace("_", "."))
    d["min_seq_sep"] = int(m.group(6))
    rows.append(d)
df = pd.concat(rows, ignore_index=True)
cols = ["ensemble", "action", "mode", "frac", "strength", "seed", "alpha", "gamma", "min_seq_sep",
        "k_pairs", "delta_mi_mean"]
df = df[cols].sort_values(cols[:9]).reset_index(drop=True)
df.to_csv(args.out, index=False)

print(f"{len(df)} alignments -> {args.out}\n")
t = (df.groupby(["ensemble", "action", "strength"]).delta_mi_mean
       .agg(mean="mean", share_positive=lambda x: (x > 0).mean()).round(3))
print(t.to_string())
