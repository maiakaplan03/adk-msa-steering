#!/usr/bin/env python3
"""
generate_phase1_sweep.py
========================

Drives the full Phase-1 MSA-editing sweep defined in HANDOFF §9 and §10.
Calls modifier_v2.modify() across every (column, theory, direction, dose)
combination and writes each resulting .a3m into the appropriate theory
folder, ready to be fed straight into the Boltz-2 prediction pipeline.

Grid:
    columns    = 50 (negative control), 64, 136, 185, 190, 202, 206
    theories   = A (amplification), B (break/shuffle), C (install modal AA)
                 — Col 50 participates only in Theory B (see HANDOFF §9).
                 — Theory A and Theory C are skipped for col 50.
    directions = closed, open
    doses      = 25, 50, 75, 100 (% of editable rows)

Total file counts (per HANDOFF §9):
    Theory A : 6 cols × 2 dir × 4 dose = 48   (col 50 skipped)
    Theory B : 7 cols × 2 dir × 4 dose = 56   (both directions confirm symmetry)
    Theory C : 6 cols × 2 dir × 4 dose = 48   (col 50 skipped)
    ----------------------------------------
                                       152 .a3m files

Output layout:
    outputs_msas/
        theory_A_amplification/
            v2__col-064__theory-A__dir-closed__dose-025.a3m
            v2__col-064__theory-A__dir-closed__dose-050.a3m
            ...
        theory_B_break/
            ...
        theory_C_install/
            ...
    manifest.csv                # one row per generated file, with all parameters

CLI:
    python generate_phase1_sweep.py                     # full sweep
    python generate_phase1_sweep.py --dry-run           # list what would run
    python generate_phase1_sweep.py --overwrite         # regenerate existing
    python generate_phase1_sweep.py --columns 64 190    # subset
    python generate_phase1_sweep.py --theories A C      # subset
    python generate_phase1_sweep.py --doses 25 100      # subset
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from typing import Dict, Iterable, List, Sequence, Tuple

# Make the modifier importable whether this script is run from this folder
# or from elsewhere.
HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import modifier_v2  # noqa: E402
from modifier_v2 import (
    DOSE_CHOICES,
    THEORY_CHOICES,
    DIRECTION_CHOICES,
    modify,
    standard_filename,
)


# ---------------------------------------------------------------------------
# Default paths — the canonical AdK inputs
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DEFAULT_ORIGINAL  = os.path.join(REPO_ROOT, "repo-folder", "AdK", "Frustration-Energy-space", "MSA.a3m")
DEFAULT_CLOSED500 = os.path.join(REPO_ROOT, "repo-folder", "AdK", "AdK-mixing", "Closed500.a3m")
DEFAULT_OPEN100   = os.path.join(REPO_ROOT, "repo-folder", "AdK", "AdK-mixing", "Open100.a3m")

OUTPUT_ROOT = os.path.join(HERE, "outputs_msas")
MANIFEST_PATH = os.path.join(HERE, "manifest.csv")

# Theory → output subfolder name. These directories already exist on disk.
THEORY_DIR = {
    "A": "theory_A_amplification",
    "B": "theory_B_break",
    "C": "theory_C_install",
}


# ---------------------------------------------------------------------------
# The Phase-1 grid (HANDOFF §9)
# ---------------------------------------------------------------------------
# (column_1based, domain_label, note)
PHASE1_COLUMNS: List[Tuple[int, str, str]] = [
    (50,  "NMP",   "negative control — Theory C skipped"),
    (64,  "CORE",  "highest-JS in JS-top set; mode V→I"),
    (136, "LID",   "Phase-0 sensitive; Theory C fallback E↔V (low freq)"),
    (185, "CORE2", "high JS; E→S switch"),
    (190, "CORE",  "highest JS overall (0.140); N→A"),
    (202, "CORE",  "Phase-0 sensitive; Theory C fallback I↔L"),
    (206, "CORE",  "high JS; F→R (aromatic↔charged, dramatic)"),
]
COLUMN_DOMAIN = {col: dom for (col, dom, _) in PHASE1_COLUMNS}
COLUMN_NOTE = {col: note for (col, _, note) in PHASE1_COLUMNS}
ALL_COLUMNS = [c for (c, _, _) in PHASE1_COLUMNS]

# HANDOFF §9 — col 50 is the negative control and only participates in Theory B.
# Theory A and Theory C are skipped for col 50.
#   "Theory A (amplification): apply to all 6 columns"        → skip col 50
#   "Theory B (break/shuffle): all 6 columns + negative control" → include col 50
#   "Theory C (install modal AA): columns 64,190,206,185,202,136 (not res 50)" → skip col 50
THEORY_SKIP_COLUMNS = {
    "A": {50},
    "B": set(),
    "C": {50},
}


# ---------------------------------------------------------------------------
# Grid expansion
# ---------------------------------------------------------------------------
def enumerate_grid(
    columns: Sequence[int],
    theories: Sequence[str],
    directions: Sequence[str],
    doses: Sequence[int],
) -> List[Dict]:
    """Expand the requested grid, applying the Theory-C skip rule for col 50."""
    plan: List[Dict] = []
    for col in columns:
        for theory in theories:
            if col in THEORY_SKIP_COLUMNS.get(theory, set()):
                continue
            for direction in directions:
                for dose in doses:
                    plan.append({
                        "column_1based": col,
                        "domain":        COLUMN_DOMAIN.get(col, "?"),
                        "theory":        theory,
                        "direction":     direction,
                        "dose_pct":      dose,
                    })
    return plan


def output_path_for(item: Dict) -> str:
    """Where this run's .a3m lives. Flat per-theory folder; filename is self-describing."""
    folder = os.path.join(OUTPUT_ROOT, THEORY_DIR[item["theory"]])
    fname = standard_filename(
        item["column_1based"], item["theory"], item["direction"], item["dose_pct"]
    )
    return os.path.join(folder, fname)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
MANIFEST_FIELDS = [
    "output_path",
    "filename",
    "theory",
    "column_1based",
    "domain",
    "direction",
    "dose_pct",
    "n_rows_total",
    "n_rows_R_j_editable",
    "n_rows_edited",
    "fraction_of_editable_rows_edited",
    "theory_C_target_aa",
    "theory_C_target_reason",
    "msa_original",
    "msa_closed500",
    "msa_open100",
    "global_seed",
    "column_note",
]


def run_sweep(
    plan: Iterable[Dict],
    msa_original: str,
    msa_closed500: str,
    msa_open100: str,
    overwrite: bool,
    dry_run: bool,
    manifest_path: str,
) -> Tuple[int, int]:
    """Execute the plan. Returns (n_generated, n_skipped_existing)."""
    plan = list(plan)
    n_total = len(plan)
    n_generated = 0
    n_skipped = 0
    manifest_rows: List[Dict] = []

    t0 = time.time()
    for idx, item in enumerate(plan, start=1):
        out = output_path_for(item)
        item["output_path"] = out
        item["filename"] = os.path.basename(out)
        item["msa_original"] = msa_original
        item["msa_closed500"] = msa_closed500
        item["msa_open100"] = msa_open100
        item["column_note"] = COLUMN_NOTE.get(item["column_1based"], "")

        already_exists = os.path.exists(out)
        if already_exists and not overwrite:
            n_skipped += 1
            print(f"  [{idx:3d}/{n_total}] SKIP (exists)  {os.path.basename(out)}")
            # Still need to populate manifest with stable params even when skipping.
            item.setdefault("n_rows_total", "")
            item.setdefault("n_rows_R_j_editable", "")
            item.setdefault("n_rows_edited", "")
            item.setdefault("fraction_of_editable_rows_edited", "")
            item.setdefault("theory_C_target_aa", "")
            item.setdefault("theory_C_target_reason", "")
            item.setdefault("global_seed", modifier_v2.GLOBAL_SEED)
            manifest_rows.append(item)
            continue

        if dry_run:
            print(f"  [{idx:3d}/{n_total}] DRY-RUN        {os.path.basename(out)}")
            continue

        summary = modify(
            msa_original=msa_original,
            msa_closed500=msa_closed500,
            msa_open100=msa_open100,
            column_1based=item["column_1based"],
            theory=item["theory"],
            direction=item["direction"],
            dose_pct=item["dose_pct"],
            output_path=out,
        )
        # Merge engine-reported facts into the manifest row.
        for k in (
            "n_rows_total",
            "n_rows_R_j_editable",
            "n_rows_edited",
            "fraction_of_editable_rows_edited",
            "global_seed",
        ):
            item[k] = summary[k]
        item["theory_C_target_aa"] = summary.get("theory_C_target_aa", "")
        item["theory_C_target_reason"] = summary.get("theory_C_target_reason", "")
        manifest_rows.append(item)
        n_generated += 1

        if idx % 25 == 0 or idx == n_total:
            elapsed = time.time() - t0
            rate = idx / elapsed if elapsed > 0 else 0.0
            print(f"  [{idx:3d}/{n_total}] WROTE          {os.path.basename(out)}   ({rate:.1f}/s)")

    if not dry_run:
        os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
        with open(manifest_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS, extrasaction="ignore")
            writer.writeheader()
            for row in manifest_rows:
                writer.writerow(row)
        print(f"\n  manifest written: {manifest_path}  ({len(manifest_rows)} rows)")

    return n_generated, n_skipped


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--msa-original",  default=DEFAULT_ORIGINAL,
                   help=f"Path to the full AdK MSA (default: {DEFAULT_ORIGINAL})")
    p.add_argument("--msa-closed500", default=DEFAULT_CLOSED500,
                   help=f"Path to Closed500 (default: {DEFAULT_CLOSED500})")
    p.add_argument("--msa-open100",   default=DEFAULT_OPEN100,
                   help=f"Path to Open100 (default: {DEFAULT_OPEN100})")
    p.add_argument("--columns", type=int, nargs="+", default=ALL_COLUMNS,
                   help=f"Columns to include (default: {ALL_COLUMNS})")
    p.add_argument("--theories", nargs="+", default=list(THEORY_CHOICES),
                   choices=list(THEORY_CHOICES),
                   help=f"Theories to include (default: {list(THEORY_CHOICES)})")
    p.add_argument("--directions", nargs="+", default=list(DIRECTION_CHOICES),
                   choices=list(DIRECTION_CHOICES),
                   help=f"Directions to include (default: {list(DIRECTION_CHOICES)})")
    p.add_argument("--doses", type=int, nargs="+", default=list(DOSE_CHOICES),
                   choices=list(DOSE_CHOICES),
                   help=f"Doses to include (default: {list(DOSE_CHOICES)})")
    p.add_argument("--overwrite", action="store_true",
                   help="Regenerate files that already exist on disk")
    p.add_argument("--dry-run", action="store_true",
                   help="List what would be written; do not write")
    p.add_argument("--manifest", default=MANIFEST_PATH,
                   help=f"Output manifest CSV path (default: {MANIFEST_PATH})")
    return p


def main(argv: List[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)

    # Validate columns
    unknown_cols = [c for c in args.columns if c not in ALL_COLUMNS]
    if unknown_cols:
        print(f"warning: column(s) {unknown_cols} are not part of the Phase-1 candidate set "
              f"{ALL_COLUMNS}; they will still be processed but check intent.",
              file=sys.stderr)

    plan = enumerate_grid(args.columns, args.theories, args.directions, args.doses)

    # Pre-flight summary
    by_theory: Dict[str, int] = {t: 0 for t in THEORY_CHOICES}
    for item in plan:
        by_theory[item["theory"]] += 1
    print("=" * 62)
    print("Phase-1 sweep plan")
    print("-" * 62)
    print(f"  columns      : {args.columns}")
    print(f"  theories     : {args.theories}")
    print(f"  directions   : {args.directions}")
    print(f"  doses        : {args.doses}")
    print(f"  total runs   : {len(plan)}")
    for t in THEORY_CHOICES:
        print(f"      theory {t}: {by_theory[t]}")
    print(f"  overwrite    : {args.overwrite}")
    print(f"  dry-run      : {args.dry_run}")
    print(f"  output root  : {OUTPUT_ROOT}")
    print(f"  manifest     : {args.manifest}")
    print("=" * 62)

    n_generated, n_skipped = run_sweep(
        plan,
        msa_original=args.msa_original,
        msa_closed500=args.msa_closed500,
        msa_open100=args.msa_open100,
        overwrite=args.overwrite,
        dry_run=args.dry_run,
        manifest_path=args.manifest,
    )

    print("-" * 62)
    print(f"  generated : {n_generated}")
    print(f"  skipped   : {n_skipped} (existing files; pass --overwrite to redo)")
    print(f"  total     : {len(plan)}")
    print("=" * 62)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
