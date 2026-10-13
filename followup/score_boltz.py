#!/usr/bin/env python3
"""
score_boltz.py
==============
Computes TM_Closed, TM_Open, NMP_Angle, LID_Angle for each Boltz .cif in
nuevos_MSA/v2_phase1/boltz_results/ and writes one tidy CSV.

Self-contained: only depends on numpy. CIF parsing is done manually; TM-score
uses iterative Kabsch superposition (Zhang & Skolnick 2004), same general
recipe used by TMscore for pre-aligned residue-corresponding chains.

References
----------
closedCA.pdb  (CA-only 1AKE, 214 residues, E. coli AdK numbering)
openCA.pdb    (CA-only 4AKE, 214 residues, E. coli AdK numbering)
"""
from __future__ import annotations
import os, re, sys, csv, time
from glob import glob
import numpy as np

from pathlib import Path

# Repo-relative paths. Override with env vars to point at data kept outside
# the repo (the 760 predicted .cif files are not distributed here).
REPO   = Path(__file__).resolve().parents[1]
INPUTS = Path(os.environ.get("ADK_INPUTS", REPO / "data" / "AFF2024"))
BOLTZ  = os.environ.get("ADK_STRUCTURES", str(REPO / "data" / "boltz_results"))
REF_CLOSED = f"{INPUTS}/AdK/Frustration-Energy-space/closedCA.pdb"
REF_OPEN   = f"{INPUTS}/AdK/Frustration-Energy-space/openCA.pdb"
OUT_CSV    = str(REPO / "followup" / "results" / "per_structure_metrics.csv")

FNAME_RE = re.compile(
    r"v2__col-(\d+)__theory-([ABC])__dir-(closed|open)__dose-(\d+)_seed(\d+)\.cif"
)

# Domain definitions: 0-based residue index range (inclusive), matching
# angle_analysis_plot_fixed.py / repo's pdb2angles.py
DOMAIN_DEFS = {
    "corelid": (114, 124),
    "core":    (89, 99),
    "nmp":     (34, 54),
    "core2":   (178, 184),
    "lid":     (124, 152),
}


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def parse_pdb_ca(path: str) -> np.ndarray:
    """Read a CA-only PDB. Returns (N, 3) array of CA coordinates in residue order."""
    coords = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            # Fixed-column PDB format
            name = line[12:16].strip()
            if name != "CA":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords.append((x, y, z))
    return np.asarray(coords, dtype=np.float64)


def parse_cif_ca(path: str) -> np.ndarray:
    """
    Read CA coordinates from a ModelCIF file (Boltz output).
    The _atom_site loop has 19 columns in order:
      group_PDB id type_symbol label_atom_id label_alt_id label_comp_id
      label_seq_id auth_seq_id ins_code label_asym_id Cartn_x Cartn_y Cartn_z
      occupancy label_entity_id auth_asym_id auth_comp_id B_iso pdbx_PDB_model_num
    """
    coords = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            atom_name = parts[3]
            if atom_name != "CA":
                continue
            try:
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
            except (ValueError, IndexError):
                continue
            coords.append((x, y, z))
    return np.asarray(coords, dtype=np.float64)


# ---------------------------------------------------------------------------
# Kabsch + TM-score
# ---------------------------------------------------------------------------
def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute optimal rotation R and translation t that brings P onto Q (P @ R + t ≈ Q).
    P, Q are (N, 3). Returns (R, t).
    """
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = U @ D @ Vt
    R = R  # 3x3
    t = Q.mean(axis=0) - P.mean(axis=0) @ R
    return R, t


def tm_score(P: np.ndarray, Q: np.ndarray, L_target: int | None = None,
             max_iter: int = 20) -> float:
    """
    TM-score with iterative refinement, for residue-aligned CA chains of equal length.

    Standard TMscore (pre-aligned) recipe:
      1. Initial Kabsch over all matched residues.
      2. Compute per-residue distances; include only residues with d < d_cut in the
         next Kabsch round. d_cut starts at d_0 and grows if too few residues fit.
      3. Iterate until set of included residues stops changing.
      4. TM = (1/L_target) * sum 1/(1 + (d_i/d_0)^2) over ALL L_target residues.

    P (model) and Q (reference) are (N, 3), N == L_target.
    """
    assert P.shape == Q.shape
    L = P.shape[0]
    if L_target is None:
        L_target = L
    # d_0 per Zhang-Skolnick 2004
    if L_target >= 30:
        d0 = 1.24 * (L_target - 15) ** (1.0 / 3.0) - 1.8
    else:
        d0 = 0.5
    d0 = max(d0, 0.5)

    # Track best score across initial-window seeds, like TMscore does.
    best_tm = -1.0
    # Use several initial alignment subsets: full chain, plus three sliding windows.
    # This helps if the prediction is locally good but globally distorted.
    L_init_set = [L]
    # Several seed windows
    for w in (L // 2, L // 3, 20):
        if w < 4 or w >= L:
            continue
        L_init_set.append(w)

    for L_init in L_init_set:
        starts = list(range(0, L - L_init + 1, max(1, L_init // 2))) or [0]
        if L_init == L:
            starts = [0]
        for s in starts:
            sel = np.arange(s, s + L_init)
            R, t = kabsch(P[sel], Q[sel])
            P_al = P @ R + t

            d_cut = d0 + 1.0
            included = sel.copy()
            for _ in range(max_iter):
                d = np.linalg.norm(P_al - Q, axis=1)
                new_inc = np.where(d < d_cut)[0]
                if len(new_inc) < 3:
                    d_cut += 0.5
                    continue
                if (len(new_inc) == len(included) and
                        np.array_equal(new_inc, included)):
                    break
                included = new_inc
                R, t = kabsch(P[included], Q[included])
                P_al = P @ R + t

            d = np.linalg.norm(P_al - Q, axis=1)
            tm = float(np.sum(1.0 / (1.0 + (d / d0) ** 2)) / L_target)
            if tm > best_tm:
                best_tm = tm
    return best_tm


# ---------------------------------------------------------------------------
# Angles (NMP-Core and LID-Core), exactly as in angle_analysis_plot_fixed.py
# ---------------------------------------------------------------------------
def domain_centroids(ca: np.ndarray) -> dict[str, np.ndarray]:
    out = {}
    for name, (lo, hi) in DOMAIN_DEFS.items():
        out[name] = ca[lo:hi + 1].mean(axis=0)
    return out


def calc_angles(ca: np.ndarray) -> tuple[float, float]:
    c = domain_centroids(ca)
    # NMP-Core: vertex at core, arms to corelid and nmp
    v1 = c["corelid"] - c["core"]
    v2 = c["nmp"]     - c["core"]
    ang_nmp = np.degrees(np.arccos(np.clip(
        np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0)))
    # LID-Core: vertex at corelid, arms to core2 and lid
    v3 = c["core2"] - c["corelid"]
    v4 = c["lid"]   - c["corelid"]
    ang_lid = np.degrees(np.arccos(np.clip(
        np.dot(v3, v4) / (np.linalg.norm(v3) * np.linalg.norm(v4)), -1.0, 1.0)))
    return float(ang_nmp), float(ang_lid)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("Loading references…")
    ref_closed = parse_pdb_ca(REF_CLOSED)
    ref_open   = parse_pdb_ca(REF_OPEN)
    assert ref_closed.shape == (214, 3), f"closedCA shape {ref_closed.shape}"
    assert ref_open.shape   == (214, 3), f"openCA   shape {ref_open.shape}"
    L = 214

    # Sanity: angles on the references themselves
    nmp_c, lid_c = calc_angles(ref_closed)
    nmp_o, lid_o = calc_angles(ref_open)
    print(f"  1AKE (closed): NMP={nmp_c:.2f}  LID={lid_c:.2f}")
    print(f"  4AKE (open)  : NMP={nmp_o:.2f}  LID={lid_o:.2f}")

    cif_files = sorted(glob(f"{BOLTZ}/theory_*/v2__*.cif"))
    print(f"Found {len(cif_files)} CIF files.")

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    rows = []
    t0 = time.time()
    bad = 0
    for i, p in enumerate(cif_files):
        fn = os.path.basename(p)
        m = FNAME_RE.match(fn)
        if not m:
            print(f"  skip (no match): {fn}")
            bad += 1
            continue
        col = int(m.group(1))
        theory = m.group(2)
        direction = m.group(3)
        dose = int(m.group(4))
        seed = int(m.group(5))

        ca = parse_cif_ca(p)
        if ca.shape != (L, 3):
            print(f"  bad shape {ca.shape} for {fn}")
            bad += 1
            continue

        tmC = tm_score(ca, ref_closed, L_target=L)
        tmO = tm_score(ca, ref_open,   L_target=L)
        nmp, lid = calc_angles(ca)

        rows.append(dict(
            filename=fn,
            column_1based=col, theory=theory, direction=direction,
            dose_pct=dose, boltz_seed=seed,
            TM_Closed=tmC, TM_Open=tmO, dTM=tmO - tmC,
            NMP_Angle=nmp, LID_Angle=lid,
        ))
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(cif_files)}  "
                  f"({(time.time()-t0):.1f}s elapsed)")

    print(f"Done in {(time.time()-t0):.1f}s.  good={len(rows)}  bad={bad}")
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {OUT_CSV}")


if __name__ == "__main__":
    main()
