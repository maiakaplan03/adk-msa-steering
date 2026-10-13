#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# check_seed_determinism.sh
#
# Goal: figure out WHY the original-MSA control produced different structures
# from the same MSA. Distinguishes three scenarios:
#   (1) diversity came from NOT fixing the seed per run (entropy),
#   (2) diversity came from one-seed-but-RNG-advancing inside a single process,
#   (3) diversity persists even with a fixed seed in separate processes
#       -> means GPU/bf16 nondeterminism is an ADDITIONAL source (handoff
#          would then be incomplete, and "subsample off" alone won't be
#          deterministic).
#
# Run on the cluster (needs the boltz env + GPU + model weights).
#
# Usage:
#   ./check_seed_determinism.sh /path/to/original_unedited.a3m
#
# Edit MSA_INPUT below or pass it as $1. The MSA must be the SAME unedited
# source MSA the control used (depth ~2651).
# ---------------------------------------------------------------------------
set -euo pipefail

MSA_INPUT="${1:?Pass the path to the original unedited MSA (.a3m or boltz yaml) as arg 1}"
SEED=123
NSUB=1024
OUT=./seed_determinism_test
mkdir -p "$OUT"

# Helper: run one prediction into a named subdir.
# Always subsample_msa=true, num_subsampled_msa=1024 to MATCH the control provenance.
run () {  # args: out_subdir  [extra boltz flags...]
  local sub="$1"; shift
  rm -rf "$OUT/$sub"; mkdir -p "$OUT/$sub"
  boltz predict "$MSA_INPUT" \
    --out_dir "$OUT/$sub" \
    --subsample_msa --num_subsampled_msa "$NSUB" \
    "$@"
}

echo "############################################################"
echo "# TEST A: SAME explicit seed, TWO SEPARATE processes"
echo "#   Expect IDENTICAL if randperm-with-fixed-seed is the only"
echo "#   stochasticity. DIFFERENT => GPU/bf16 nondeterminism (scenario 3)."
echo "############################################################"
run A_seed_run1 --seed "$SEED"
run A_seed_run2 --seed "$SEED"

echo "############################################################"
echo "# TEST B: NO seed, TWO SEPARATE processes"
echo "#   Expect DIFFERENT (entropy seeding) -> scenario 1 reproduced."
echo "############################################################"
run B_noseed_run1
run B_noseed_run2

echo
echo "=== Comparing CA coordinates (definitive) ==="
python3 - "$OUT" <<'PY'
import sys, glob, os
import numpy as np
try:
    import gemmi
except ImportError:
    sys.exit("Install gemmi in this env to compare structures: pip install gemmi")

root = sys.argv[1]
def ca_coords(cif):
    st = gemmi.read_structure(cif); st.setup_entities()
    xs = []
    for model in st:
        for chain in model:
            for res in chain:
                a = res.find_atom("CA", "*")
                if a: xs.append([a.pos.x, a.pos.y, a.pos.z])
        break
    return np.array(xs)

def find_cif(subdir):
    hits = glob.glob(os.path.join(root, subdir, "**", "*_model_0.cif"), recursive=True)
    return hits[0] if hits else None

def cmp(label, s1, s2):
    c1, c2 = find_cif(s1), find_cif(s2)
    if not c1 or not c2:
        print(f"{label}: missing output ({c1}, {c2})"); return
    a, b = ca_coords(c1), ca_coords(c2)
    if a.shape != b.shape:
        print(f"{label}: shape mismatch {a.shape} vs {b.shape}"); return
    rmsd = float(np.sqrt(((a-b)**2).sum(axis=1).mean()))
    maxd = float(np.sqrt(((a-b)**2).sum(axis=1)).max())
    verdict = "IDENTICAL" if maxd < 1e-3 else ("near-identical" if rmsd < 0.1 else "DIFFERENT")
    print(f"{label}: CA-RMSD={rmsd:.4f} A, max atom dev={maxd:.4f} A  -> {verdict}")

print()
cmp("TEST A (same seed, 2 processes) ", "A_seed_run1", "A_seed_run2")
cmp("TEST B (no seed, 2 processes)   ", "B_noseed_run1", "B_noseed_run2")
print()
print("Interpretation:")
print("  A IDENTICAL + B DIFFERENT  -> diversity was lack-of-fixed-seed (scenarios 1/2).")
print("                                Handoff complete; 'subsample off' OR per-run seed = deterministic.")
print("  A DIFFERENT                -> GPU/bf16 nondeterminism is a real extra source (scenario 3).")
print("                                Even fixed seeds won't be deterministic; need deterministic kernels.")
PY
