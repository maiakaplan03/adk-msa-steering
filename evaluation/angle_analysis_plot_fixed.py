#!/usr/bin/env python3
"""
Patched version of angle_analysis_plot.py.

Fixes vs. the original:
  1. LID-Core angle is now computed with the vertex at `corelid` (arms to
     `core2` and `lid`), matching the repo's pdb2angles.py / traj2angle.py.
     The original placed the vertex at `core2`, which is a different angle
     and was the reason the y-axis range disagreed with the paper.
  2. 1AKE.pdb / 4AKE.pdb are detected as reference structures and drawn
     as distinct markers (red filled circle for 1ake, red open circle
     for 4ake), like the paper.
  3. Default --dir now resolves relative to *this script's location*,
     so running from any cwd works as long as `structures/` sits next
     to the script.
  4. Renamed the shadowed `c_lid` -> `c_corelid` to avoid confusion with
     the actual `lid` centroid.
"""

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mdtraj as md

# Domain definitions (0-based residue indices for E. coli AdK).
# Same as the repo's pdb2angles.py.
DOMAIN_DEFS = {
    "corelid": (114, 124),
    "core":    (89, 99),
    "nmp":     (34, 54),
    "core2":   (178, 184),
    "lid":     (124, 152),
}

# Filenames (lowercased) plotted as references rather than beta-colored points.
REFERENCE_FILES = {"1ake.pdb", "4ake.pdb"}


def calculate_angles(traj):
    """
    Calculates NMP-Core and LID-Core angles using the same geometric
    definitions as the repo:
      - NMP-Core: vertex at `core`,    arms to `corelid` and `nmp`
      - LID-Core: vertex at `corelid`, arms to `core2`   and `lid`
    Returns (ang_nmp, ang_lid) in degrees, or (None, None) on failure.
    """
    centroids = {}
    ca_atoms = traj.topology.select("name CA")
    for name, (lo, hi) in DOMAIN_DEFS.items():
        residue_atoms = traj.topology.select(f"resid >= {lo} and resid <= {hi}")
        indices = np.intersect1d(ca_atoms, residue_atoms)
        if len(indices) == 0:
            return None, None
        centroids[name] = np.mean(traj.xyz[:, indices, :], axis=1)

    c_corelid = centroids["corelid"][0]
    c_core    = centroids["core"][0]
    c_core2   = centroids["core2"][0]
    c_nmp     = centroids["nmp"][0]
    c_lid     = centroids["lid"][0]

    # NMP-Core: angle at `core`, between corelid and nmp
    v1 = c_corelid - c_core
    v2 = c_nmp     - c_core
    cos_nmp = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    ang_nmp = np.degrees(np.arccos(np.clip(cos_nmp, -1.0, 1.0)))

    # LID-Core: angle at `corelid`, between core2 and lid  (FIXED)
    v3 = c_core2 - c_corelid
    v4 = c_lid   - c_corelid
    cos_lid = np.dot(v3, v4) / (np.linalg.norm(v3) * np.linalg.norm(v4))
    ang_lid = np.degrees(np.arccos(np.clip(cos_lid, -1.0, 1.0)))

    return ang_nmp, ang_lid


def parse_beta(fname):
    """
    Parses beta value from filenames like 'beta_pos0p51_...' or
    'beta_neg0p38_...'. Returns 0.0 if the pattern is not present.
    """
    match = re.search(r"beta_(pos|neg)(\d+)p(\d+)", fname)
    if match:
        sign = 1.0 if match.group(1) == "pos" else -1.0
        val = float(f"{match.group(2)}.{match.group(3)}")
        return sign * val
    return 0.0


def main():
    # Resolve defaults relative to *this script's location*, so the script
    # works no matter what cwd you launch it from, as long as the
    # `structures` folder sits next to the script.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_struct_dir = os.path.join(script_dir, "structures")
    default_out = os.path.join(script_dir, "beta_gradient_analysis.png")

    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default=default_struct_dir,
                        help="Path to PDB folder (default: <script_dir>/structures)")
    parser.add_argument("--out", default=default_out,
                        help="Output PNG path (default: <script_dir>/beta_gradient_analysis.png)")
    args = parser.parse_args()

    struct_dir = os.path.abspath(args.dir)
    if not os.path.isdir(struct_dir):
        raise SystemExit(f"Structure folder not found: {struct_dir}")

    results = []
    references = {}  # filename(lower) -> (nmp, lid)

    pdb_files = [f for f in os.listdir(struct_dir) if f.endswith(".pdb")]
    print(f"Reading from: {struct_dir}")
    print(f"Found {len(pdb_files)} structures. Processing...")

    for fname in pdb_files:
        try:
            traj = md.load(os.path.join(struct_dir, fname))
            nmp, lid = calculate_angles(traj)
            if nmp is None:
                continue

            if fname.lower() in REFERENCE_FILES:
                references[fname.lower()] = (nmp, lid)
            else:
                results.append({
                    "fname": fname,
                    "beta":  parse_beta(fname),
                    "nmp":   nmp,
                    "lid":   lid,
                })
        except Exception as e:
            print(f"Error processing {fname}: {e}")

    if not results and not references:
        print("No valid data points calculated.")
        return

    df = pd.DataFrame(results)

    # --- Plot ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

    if not df.empty:
        # Center colormap around 0 so red = positive beta, blue = negative.
        vmax = df["beta"].abs().max() if df["beta"].abs().max() > 0 else 1.0
        vmin = -vmax

        sc = ax.scatter(
            df["nmp"], df["lid"],
            c=df["beta"],
            cmap="coolwarm",
            vmin=vmin, vmax=vmax,
            s=40, alpha=0.7, edgecolors="none",
        )
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label("Beta Value (Gradient)")

    # Overlay 1ake / 4ake as references, mimicking the paper's style.
    if "1ake.pdb" in references:
        nmp, lid = references["1ake.pdb"]
        ax.scatter(nmp, lid, s=180, c="red", marker="o",
                   edgecolors="darkred", linewidths=1.5,
                   zorder=5, label="1ake (closed)")
        ax.annotate("1ake", (nmp, lid), xytext=(6, -14),
                    textcoords="offset points", color="red",
                    fontsize=12, fontweight="bold")

    if "4ake.pdb" in references:
        nmp, lid = references["4ake.pdb"]
        ax.scatter(nmp, lid, s=180, facecolors="none", edgecolors="red",
                   linewidths=2.0, marker="o",
                   zorder=5, label="4ake (open)")
        ax.annotate("4ake", (nmp, lid), xytext=(6, 6),
                    textcoords="offset points", color="red",
                    fontsize=12, fontweight="bold")

    ax.set_xlabel("Angle NMP-Core (\u00b0)")
    ax.set_ylabel("Angle LID-Core (\u00b0)")
    ax.set_title(f"Conformational Analysis (n={len(results)})")

    if references:
        ax.legend(loc="best", framealpha=0.9)

    plt.tight_layout()
    plt.savefig(args.out)
    plt.show()
    print(f"Saved plot to {args.out}")


if __name__ == "__main__":
    main()
