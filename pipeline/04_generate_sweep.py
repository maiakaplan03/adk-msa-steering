#!/usr/bin/env python3
"""
04_generate_sweep.py — "4- generamos un barrido de MSAs artificiales de AdK".

Maia's sweep loop (notebook "BLOQUE 2 — SWEEP / BARRIDO") as a command-line script.
Every combination of the grid calls msa_editor.make_progressive_msas_weighted() and the
resulting alignments are copied into one master folder with the parameters encoded in the
file name, e.g.

  msa_degree_frac0.30_str0.30_seed11__ACT-promote_joint__MODE-degree__A-0_05__G-2__MINSEP-12.a3m

Two grids are provided:

  --preset ensembles     (default) the grid behind the four 324-alignment ensembles that
                         were folded with Boltz-2 and scored (results/sweep_evaluation):
                         frac 0.3/0.6/0.9, strength 0.3/0.6/0.9, seeds 11/22/33,
                         alpha 0.05/0.1, gamma 2/3/4, min_seq_sep 8/12, mode degree.
                         3*3*3*2*3*2 = 324 alignments per (contact set, action).

  --preset presentation  the grid shown in the TP presentation (slide 4): actions promote
                         and promote_joint, frac 0.2/0.4/0.6/0.8, strength 0.4/0.6/0.8,
                         same seeds/alpha/gamma/min_seq_sep. 864 alignments. File names
                         carry the extra "__DF-0" tag (distance filter off) that the
                         presentation-era files had.

Ensembles (contact set x action), named as the original output folders:

  A_minus_B                   contacts lost from 1AKE to 4AKE (closed-specific, 66)  promote_joint
  B_minus_A                   contacts gained from 1AKE to 4AKE (open-specific, 57)  promote_joint
  A_minus_B_control_negativo  closed-specific contacts                               break
  B_minus_A_control_negativo  open-specific contacts                                 break

Examples
  python pipeline/04_generate_sweep.py --msa data/AFF2024/AdK/Frustration-Energy-space/MSA.a3m \
         --ensemble all --out-root data/sweep_msas
  shasum -a 256 -c data/sweep_msa_checksums.sha256      # run from data/sweep_msas

The defaults reproduce the 1,296 alignments that were run, byte for byte.
"""
import argparse
import itertools
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from msa_editor import make_progressive_msas_weighted, copy_a3m_tree_to_master, _fmt  # noqa: E402

CONTACTS = os.path.join(REPO, "results", "contacts")
ENSEMBLES = {
    "A_minus_B":                  ("delta_contacts_lost_A_minus_B.csv",   "promote_joint"),
    "B_minus_A":                  ("delta_contacts_gained_B_minus_A.csv", "promote_joint"),
    "A_minus_B_control_negativo": ("delta_contacts_lost_A_minus_B.csv",   "break"),
    "B_minus_A_control_negativo": ("delta_contacts_gained_B_minus_A.csv", "break"),
}

PRESETS = {
    "ensembles": dict(
        actions=None,  # taken from the ensemble
        modes=["degree"],
        pair_fracs=[0.3, 0.6, 0.9],
        strengths=[0.3, 0.6, 0.9],
        seeds=[11, 22, 33],
        alphas=[0.05, 0.1],
        gammas=[2.0, 3.0, 4.0],
        min_seq_seps=[8, 12],
        df_tag=False,
    ),
    "presentation": dict(
        actions=["promote", "promote_joint"],
        modes=["degree"],
        pair_fracs=[0.2, 0.4, 0.6, 0.8],
        strengths=[0.4, 0.6, 0.8],
        seeds=[11, 22, 33],
        alphas=[0.05, 0.1],
        gammas=[2.0, 3.0, 4.0],
        min_seq_seps=[8, 12],
        df_tag=True,
    ),
}


def run_sweep(msa, contacts_csv, master_dir, grid, actions, overlap="sequential", keep_runs=False):
    base_out_dir = master_dir + "__runs"
    os.makedirs(base_out_dir, exist_ok=True)
    os.makedirs(master_dir, exist_ok=True)
    total = 0
    for action, mode, alpha, gamma, min_sep in itertools.product(
            actions, grid["modes"], grid["alphas"], grid["gammas"], grid["min_seq_seps"]):
        run_out_dir = os.path.join(
            base_out_dir,
            f"run__act-{action}__mode-{mode}__a-{_fmt(alpha)}__g-{_fmt(gamma)}__minsep-{min_sep}")
        make_progressive_msas_weighted(
            msa_in_path=msa, contacts_csv=contacts_csv, base_out_dir=run_out_dir,
            mode=mode, pair_fracs=grid["pair_fracs"], strength_grid=grid["strengths"],
            seeds=grid["seeds"], action=action, ref_index=0, one_based=True,
            alpha=alpha, power_gamma=gamma, min_seq_sep=min_sep, overlap=overlap,
            use_dist_filter=False, max_dist_for_gain=None)
        suffix = (f"__ACT-{action}__MODE-{mode}__A-{_fmt(alpha)}__G-{_fmt(gamma)}"
                  f"__MINSEP-{min_sep}" + ("__DF-0" if grid["df_tag"] else ""))
        copy_a3m_tree_to_master(src_root=run_out_dir, dst_root=master_dir, extra_suffix=suffix)
        total += 1
        print(f"  run {total}: {suffix}")
    if not keep_runs:
        shutil.rmtree(base_out_dir, ignore_errors=True)
    n = len([f for f in os.listdir(master_dir) if f.endswith(".a3m")])
    print(f"  -> {n} alignments in {master_dir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--msa", required=True, help="source alignment (Guan et al. MSA.a3m, depth 2651)")
    ap.add_argument("--ensemble", default="all", choices=list(ENSEMBLES) + ["all"])
    ap.add_argument("--contacts-dir", default=CONTACTS, help="folder with the delta_contacts_*.csv files")
    ap.add_argument("--out-root", default=os.path.join(REPO, "data", "sweep_msas"))
    ap.add_argument("--preset", default="ensembles", choices=list(PRESETS))
    ap.add_argument("--overlap", default="sequential", choices=["sequential", "lock"],
                    help="'sequential' = original behaviour (reproduces the checksums); "
                         "'lock' = review fix, one write per column (see docs/review_notes.md)")
    ap.add_argument("--keep-runs", action="store_true", help="keep per-run folders with summary/report CSVs")
    args = ap.parse_args()

    grid = PRESETS[args.preset]
    names = list(ENSEMBLES) if args.ensemble == "all" else [args.ensemble]
    for name in names:
        csv_name, action = ENSEMBLES[name]
        actions = grid["actions"] or [action]
        if grid["actions"] and name.endswith("control_negativo"):
            actions = ["break"]
        master = os.path.join(args.out_root, f"all_runs_a3m_{name}")
        print(f"[{name}] contacts={csv_name} actions={actions} preset={args.preset} overlap={args.overlap}")
        run_sweep(args.msa, os.path.join(args.contacts_dir, csv_name), master, grid, actions,
                  overlap=args.overlap, keep_runs=args.keep_runs)


if __name__ == "__main__":
    main()
