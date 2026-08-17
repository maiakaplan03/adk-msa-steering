#!/usr/bin/env python3
"""
modifier_v2.py
==============

Deterministic, single-column MSA modifier for the Phase-1 AdK steering
experiments. Implements the substitution recipe locked in HANDOFF §9.

Three independent "theories" — each tests a different mechanism by which
the predictor might be reading the MSA at a given column:

    Theory A (amplification)
        For each selected row, sample a new amino acid at column j from
        the target subfamily's column distribution (Closed500 or Open100).
        Changes the per-column amino acid COMPOSITION toward the target.

    Theory B (anchor removal / break)
        For each selected row, sample a new amino acid at column j from
        the multiset of residues currently in that column (V). Disrupts
        the row-pairing (and thus MI / co-evolution) without changing
        the per-column amino acid composition in expectation.

    Theory C (switch installation)
        Install a single deterministic amino acid (modal AA of the target
        subfamily, with 2nd-modal fallback when modes agree) at column j
        in every selected row.

Invariants (HANDOFF §11) enforced at every call:
    * Sequence length is exactly 214 columns (the AdK MSAs).
    * The query row (row 0) is NEVER edited.
    * Gap positions ('-' or lowercase a3m inserts) are NEVER introduced
      or removed.
    * Row selection is deterministic with seed = 0, and doses 25/50/75/100
      are nested by construction (rows25 ⊂ rows50 ⊂ rows75 ⊂ rows100).
    * Per-row RNG seeds (Theories A and B) make the residue actually
      assigned to a given row stable across dose tiers as well — so the
      lower-dose MSAs are literal subsets of the higher-dose MSAs at the
      level of individual substitutions, not just row selections.

Library usage:
    from modifier_v2 import modify, standard_filename
    summary = modify(msa_original, msa_closed500, msa_open100,
                     column_1based=64, theory="A",
                     direction="closed", dose_pct=25,
                     output_path="...")

CLI usage:
    python modifier_v2.py \\
        --msa-original   .../MSA.a3m \\
        --msa-closed500  .../Closed500.a3m \\
        --msa-open100    .../Open100.a3m \\
        --column 64                 # 1-based residue number
        --theory A                  # A | B | C
        --direction closed          # closed | open
        --dose 25                   # 25 | 50 | 75 | 100
        --output /path/to/out.a3m
"""

from __future__ import annotations

import argparse
import hashlib
import os
import random
import re
from collections import Counter
from typing import Dict, List, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Constants — drawn from HANDOFF §9 / §11
# ---------------------------------------------------------------------------
GAP_CHARS: set = {"-", "."}
LOWER_RE = re.compile(r"[a-z]")
N_COLS_EXPECTED: int = 214      # all three AdK MSAs share this column count
QUERY_ROW: int = 0              # row 0 is always the E. coli AdK query
GLOBAL_SEED: int = 0            # locked seed for row-selection shuffle

THEORY_CHOICES: Tuple[str, ...] = ("A", "B", "C")
DIRECTION_CHOICES: Tuple[str, ...] = ("closed", "open")
DOSE_CHOICES: Tuple[int, ...] = (25, 50, 75, 100)


# ---------------------------------------------------------------------------
# a3m I/O
# ---------------------------------------------------------------------------
def is_gap(ch: str) -> bool:
    """Treat '-', '.' and any lowercase letter (a3m insert) as a gap."""
    return ch in GAP_CHARS or (ch.isalpha() and ch.islower())


def read_a3m(path: str) -> Tuple[List[str], List[str]]:
    """Read a FASTA / a3m file. Returns (headers, sequences) preserving order."""
    headers: List[str] = []
    seqs: List[str] = []
    with open(path) as f:
        cur_header: str | None = None
        cur_parts: List[str] = []
        for line in f:
            line = line.rstrip("\n").rstrip("\r").strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_header is not None:
                    seqs.append("".join(cur_parts))
                cur_header = line[1:]
                headers.append(cur_header)
                cur_parts = []
            else:
                cur_parts.append(line)
        if cur_header is not None:
            seqs.append("".join(cur_parts))
    return headers, seqs


def normalize_a3m(seqs: List[str]) -> List[str]:
    """Map lowercase characters (a3m insert columns) to '-'.

    The AdK MSAs already have no inserts, but this is defensive: it puts
    every sequence on the same aligned column grid before column-level
    operations. Idempotent on already-clean input.
    """
    return [LOWER_RE.sub("-", s) for s in seqs]


def write_a3m(path: str, headers: List[str], seqs: List[str]) -> None:
    """Write a one-line-per-sequence a3m file."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as f:
        for h, s in zip(headers, seqs):
            f.write(f">{h}\n{s}\n")


def get_column(seqs: List[str], col_idx: int) -> List[str]:
    """Return the (0-based) column at col_idx across every row."""
    return [s[col_idx] for s in seqs]


def nongap_aa_counts(column: List[str]) -> Counter:
    """Counter of non-gap residues in a column."""
    return Counter(ch for ch in column if not is_gap(ch))


# ---------------------------------------------------------------------------
# Target distributions
# ---------------------------------------------------------------------------
def target_distribution(
    closed500_seqs: List[str],
    open100_seqs: List[str],
    col_idx: int,
    direction: str,
) -> Tuple[List[str], np.ndarray]:
    """Empirical AA distribution at column ``col_idx`` for the target subfamily.

    Returns (aas_sorted, probs_normalized). Used by Theory A.
    """
    pool = closed500_seqs if direction == "closed" else open100_seqs
    counts = nongap_aa_counts(get_column(pool, col_idx))
    if not counts:
        raise ValueError(
            f"Empty non-gap AA distribution at column {col_idx + 1} "
            f"in the {direction!r} subfamily."
        )
    aas = sorted(counts.keys())
    probs = np.array([counts[a] for a in aas], dtype=float)
    probs /= probs.sum()
    return aas, probs


def modal_target_aa(
    closed500_seqs: List[str],
    open100_seqs: List[str],
    col_idx: int,
    direction: str,
) -> Tuple[str, str]:
    """HANDOFF §9 Theory C target rule.

    If the closed and open modes differ → use the target's mode.
    If they agree                       → fall back to the target's 2nd-most-frequent AA.

    Returns (target_aa, reason) where ``reason`` is "mode_differs" or "fallback_2nd_mode".
    """
    c_counts = nongap_aa_counts(get_column(closed500_seqs, col_idx))
    o_counts = nongap_aa_counts(get_column(open100_seqs, col_idx))
    if not c_counts or not o_counts:
        raise ValueError(
            f"Empty AA distribution at column {col_idx + 1} in one of the subfamilies "
            f"(closed={len(c_counts)} aas, open={len(o_counts)} aas)."
        )

    m_c = c_counts.most_common(1)[0][0]
    m_o = o_counts.most_common(1)[0][0]

    if m_c != m_o:
        return (m_c if direction == "closed" else m_o), "mode_differs"

    # modes agree → use 2nd-most-frequent of the target distribution
    target_counts = c_counts if direction == "closed" else o_counts
    ranked = target_counts.most_common()
    if len(ranked) < 2:
        raise ValueError(
            f"Theory C fallback failed at column {col_idx + 1}: "
            f"target subfamily {direction!r} has only one non-gap amino acid."
        )
    return ranked[1][0], "fallback_2nd_mode"


# ---------------------------------------------------------------------------
# Row selection (HANDOFF §9)
# ---------------------------------------------------------------------------
def editable_rows(msa: List[str], col_idx: int) -> List[int]:
    """Indices i where MSA[i][col_idx] is NOT a gap AND i != QUERY_ROW.

    Editable rows at column j are R_j in the HANDOFF spec.
    """
    return [
        i for i in range(len(msa))
        if i != QUERY_ROW and not is_gap(msa[i][col_idx])
    ]


def selected_rows_for_dose(
    R_j: List[int],
    dose_pct: int,
    seed: int = GLOBAL_SEED,
) -> List[int]:
    """Shuffle R_j with `seed`, return the first ⌊dose · |R_j| / 100⌋ rows.

    Because the shuffle is deterministic and we take a prefix, doses
    25/50/75/100 are nested:  S_25 ⊂ S_50 ⊂ S_75 ⊂ S_100 = R_j.
    """
    rng = random.Random(seed)
    order = list(R_j)
    rng.shuffle(order)
    k = (dose_pct * len(order)) // 100   # explicit floor division
    return order[:k]


# ---------------------------------------------------------------------------
# Per-row deterministic RNG (for Theories A and B)
# ---------------------------------------------------------------------------
def per_row_rng(
    column_1based: int,
    row_idx: int,
    theory: str,
    direction: str,
) -> np.random.Generator:
    """A numpy Generator keyed reproducibly on (column, row, theory, direction).

    This is what makes the per-row substitutions for Theories A and B
    stable across dose tiers: the choice made for row i depends only on
    its identity, not on the size of the selection it ended up in.
    """
    key = (
        f"modifier_v2|col={column_1based}|row={row_idx}"
        f"|theory={theory}|direction={direction}"
    ).encode("utf-8")
    digest = hashlib.sha256(key).digest()
    seed = int.from_bytes(digest[:8], "big")
    return np.random.default_rng(seed)


# ---------------------------------------------------------------------------
# Theory implementations
# ---------------------------------------------------------------------------
def _replace_char(seq: str, col_idx: int, new_aa: str) -> str:
    """Return a copy of seq with seq[col_idx] replaced by new_aa."""
    return seq[:col_idx] + new_aa + seq[col_idx + 1:]


def apply_theory_A(
    seqs: List[str],
    col_idx: int,
    rows_S: List[int],
    closed500: List[str],
    open100: List[str],
    direction: str,
    column_1based: int,
) -> List[str]:
    aas, probs = target_distribution(closed500, open100, col_idx, direction)
    out = list(seqs)
    for i in rows_S:
        rng = per_row_rng(column_1based, i, "A", direction)
        new_aa = str(rng.choice(aas, p=probs))
        out[i] = _replace_char(out[i], col_idx, new_aa)
    return out


def apply_theory_B(
    seqs: List[str],
    col_idx: int,
    rows_S: List[int],
    R_j: List[int],
    direction: str,
    column_1based: int,
) -> List[str]:
    """V = multiset of residues currently at col j across the FULL editable
    pool R_j (not just the dose-selected subset). Each selected row picks
    one residue from V using its own seeded RNG.

    The `direction` argument is included only so we can label the output
    file; per HANDOFF §9 it has no actual effect on Theory B output — the
    two directions are run to confirm the symmetry.
    """
    V: List[str] = [seqs[i][col_idx] for i in R_j]
    if not V:
        raise ValueError(
            f"Theory B has empty residue pool V at column {col_idx + 1}."
        )
    out = list(seqs)
    for i in rows_S:
        rng = per_row_rng(column_1based, i, "B", direction)
        new_aa = V[int(rng.integers(0, len(V)))]
        out[i] = _replace_char(out[i], col_idx, new_aa)
    return out


def apply_theory_C(
    seqs: List[str],
    col_idx: int,
    rows_S: List[int],
    closed500: List[str],
    open100: List[str],
    direction: str,
) -> Tuple[List[str], str, str]:
    target_aa, reason = modal_target_aa(closed500, open100, col_idx, direction)
    out = list(seqs)
    for i in rows_S:
        out[i] = _replace_char(out[i], col_idx, target_aa)
    return out, target_aa, reason


# ---------------------------------------------------------------------------
# Top-level orchestration
# ---------------------------------------------------------------------------
def modify(
    msa_original: str,
    msa_closed500: str,
    msa_open100: str,
    column_1based: int,
    theory: str,
    direction: str,
    dose_pct: int,
    output_path: str,
) -> Dict:
    """Produce one edited MSA per (column, theory, direction, dose) and write it.

    Returns a summary dict including row counts and (for Theory C) the
    installed AA.
    """
    if theory not in THEORY_CHOICES:
        raise ValueError(f"theory must be one of {THEORY_CHOICES}, got {theory!r}")
    if direction not in DIRECTION_CHOICES:
        raise ValueError(f"direction must be one of {DIRECTION_CHOICES}, got {direction!r}")
    if dose_pct not in DOSE_CHOICES:
        raise ValueError(f"dose must be one of {DOSE_CHOICES}, got {dose_pct!r}")

    headers, seqs_raw = read_a3m(msa_original)
    seqs = normalize_a3m(seqs_raw)
    _, c500_raw = read_a3m(msa_closed500)
    _, o100_raw = read_a3m(msa_open100)
    c500 = normalize_a3m(c500_raw)
    o100 = normalize_a3m(o100_raw)

    # --- invariants -------------------------------------------------------
    L = len(seqs[0])
    if L != N_COLS_EXPECTED:
        raise ValueError(
            f"Expected {N_COLS_EXPECTED} columns in original MSA, got {L}"
        )
    if any(len(s) != L for s in seqs):
        raise ValueError("Ragged MSA: not all sequences in original have the same length.")
    if any(len(s) != L for s in c500):
        raise ValueError("Closed500 MSA has a different column count than original.")
    if any(len(s) != L for s in o100):
        raise ValueError("Open100 MSA has a different column count than original.")
    if not (1 <= column_1based <= L):
        raise ValueError(f"--column {column_1based} out of range [1, {L}]")

    col_idx = column_1based - 1
    R_j = editable_rows(seqs, col_idx)
    rows_S = selected_rows_for_dose(R_j, dose_pct, seed=GLOBAL_SEED)

    target_aa: str | None = None
    target_reason: str | None = None

    if theory == "A":
        new_seqs = apply_theory_A(seqs, col_idx, rows_S, c500, o100, direction, column_1based)
    elif theory == "B":
        new_seqs = apply_theory_B(seqs, col_idx, rows_S, R_j, direction, column_1based)
    elif theory == "C":
        new_seqs, target_aa, target_reason = apply_theory_C(
            seqs, col_idx, rows_S, c500, o100, direction
        )
    else:  # pragma: no cover — guarded above
        raise AssertionError(f"unreachable theory={theory}")

    # --- post-condition checks -------------------------------------------
    if new_seqs[QUERY_ROW] != seqs[QUERY_ROW]:
        raise AssertionError("Query row (row 0) was modified — invariant violated.")

    gaps_before = [is_gap(c) for c in get_column(seqs, col_idx)]
    gaps_after = [is_gap(c) for c in get_column(new_seqs, col_idx)]
    if gaps_before != gaps_after:
        raise AssertionError(
            f"Gap pattern at column {column_1based} changed — invariant violated."
        )

    edited_set = set(rows_S)
    for i in range(len(seqs)):
        if i in edited_set:
            # the only character that may differ is column col_idx
            for j in range(L):
                if j == col_idx:
                    continue
                if new_seqs[i][j] != seqs[i][j]:
                    raise AssertionError(
                        f"Row {i} was modified outside the target column "
                        f"(at column {j + 1})."
                    )
        else:
            if new_seqs[i] != seqs[i]:
                raise AssertionError(
                    f"Row {i} was modified but is not in the dose selection."
                )

    write_a3m(output_path, headers, new_seqs)

    summary: Dict = {
        "output_path": output_path,
        "msa_original": msa_original,
        "msa_closed500": msa_closed500,
        "msa_open100": msa_open100,
        "column_1based": column_1based,
        "theory": theory,
        "direction": direction,
        "dose_pct": dose_pct,
        "n_rows_total": len(seqs),
        "n_rows_R_j_editable": len(R_j),
        "n_rows_edited": len(rows_S),
        "fraction_of_editable_rows_edited": (
            len(rows_S) / len(R_j) if R_j else 0.0
        ),
        "global_seed": GLOBAL_SEED,
    }
    if theory == "C":
        summary["theory_C_target_aa"] = target_aa
        summary["theory_C_target_reason"] = target_reason
    return summary


def standard_filename(
    column_1based: int,
    theory: str,
    direction: str,
    dose_pct: int,
) -> str:
    """The HANDOFF §9 file-naming convention.

    Example: ``v2__col-064__theory-A__dir-closed__dose-025.a3m``
    """
    return (
        f"v2__col-{column_1based:03d}"
        f"__theory-{theory}"
        f"__dir-{direction}"
        f"__dose-{dose_pct:03d}.a3m"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--msa-original", required=True,
                   help="Path to the full AdK MSA (2,651 × 214). e.g. MSA.a3m")
    p.add_argument("--msa-closed500", required=True,
                   help="Path to the closed-biased subfamily MSA (501 × 214). e.g. Closed500.a3m")
    p.add_argument("--msa-open100", required=True,
                   help="Path to the open-biased subfamily MSA (101 × 214). e.g. Open100.a3m")
    p.add_argument("--column", type=int, required=True,
                   help="1-based residue number to edit (1..214)")
    p.add_argument("--theory", choices=THEORY_CHOICES, required=True,
                   help="Substitution theory to apply")
    p.add_argument("--direction", choices=DIRECTION_CHOICES, required=True,
                   help="Target subfamily for amino-acid choice")
    p.add_argument("--dose", type=int, choices=DOSE_CHOICES, required=True,
                   help="Fraction of editable rows to modify (25, 50, 75, or 100)")
    p.add_argument("--output", required=True,
                   help="Output .a3m path (parent directories created automatically)")
    return p


def main(argv: List[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    summary = modify(
        msa_original=args.msa_original,
        msa_closed500=args.msa_closed500,
        msa_open100=args.msa_open100,
        column_1based=args.column,
        theory=args.theory,
        direction=args.direction,
        dose_pct=args.dose,
        output_path=args.output,
    )
    width = max(len(k) for k in summary)
    for k, v in summary.items():
        print(f"  {k:<{width}}  {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
