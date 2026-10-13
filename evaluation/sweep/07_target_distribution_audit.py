"""
07_target_distribution_audit.py

Sanity checks before designing the deterministic modifier:
  1) Do Closed500 / Open100 share the original MSA's query and column count?
  2) At each candidate column (136, 202, 172, 155, 192, 109), how different
     is the residue distribution between Closed500 and Open100?
     Report:
        - top-3 residues + frequencies for original / Closed500 / Open100
        - Jensen-Shannon divergence between the two target distributions
        - the "Theory C" mode-residue for each target
"""

import os, numpy as np
from collections import Counter

# --- repo-relative paths -------------------------------------------------
from pathlib import Path as _Path
import os as _os
REPO   = _Path(__file__).resolve().parents[2]
RUNS   = _Path(_os.environ.get("ADK_RUNS",   REPO.parent / "nuevos_MSA" / "fuente_MSA_original" / "results"))
INPUTS = _Path(_os.environ.get("ADK_INPUTS", REPO / "data" / "AFF2024"))
# -------------------------------------------------------------------------


ORIG_PATH = str(INPUTS / "AdK" / "Frustration-Energy-space" / "MSA.a3m")
CLO_PATH  = str(INPUTS / "AdK" / "AdK-mixing" / "Closed500.a3m")
OPN_PATH  = str(INPUTS / "AdK" / "AdK-mixing" / "Open100.a3m")

def load_a3m(path):
    """Returns list of (header, sequence). Sequences kept exactly as in file."""
    out = []
    with open(path) as fh:
        h, s = None, []
        for line in fh:
            line = line.rstrip()
            if not line: continue
            if line.startswith(">"):
                if h is not None:
                    out.append((h, "".join(s)))
                h = line[1:]; s = []
            else:
                s.append(line)
        if h is not None:
            out.append((h, "".join(s)))
    return out

def strip_a3m_inserts(seq):
    """A3M: keep uppercase + '-' ; lowercase = insertion vs query, strip them."""
    return "".join(c for c in seq if not c.islower())

def col_dist(seqs, col, include_gaps=False):
    """Counter of residues at a given column."""
    c = Counter()
    for s in seqs:
        ch = s[col]
        if ch == "-" and not include_gaps:
            continue
        c[ch] += 1
    return c

def freqs(counter):
    n = sum(counter.values())
    if n == 0: return {}
    return {k: v/n for k, v in counter.items()}

def js_divergence(p, q):
    """Jensen-Shannon divergence (log base 2). 0 = identical, 1 = max."""
    keys = set(p) | set(q)
    pv = np.array([p.get(k,0) for k in keys])
    qv = np.array([q.get(k,0) for k in keys])
    m = 0.5 * (pv + qv)
    def kl(a, b):
        mask = (a > 0) & (b > 0)
        return np.sum(a[mask] * np.log2(a[mask] / b[mask]))
    return 0.5 * kl(pv, m) + 0.5 * kl(qv, m)

# --- load ---
orig_rec = load_a3m(ORIG_PATH)
clo_rec  = load_a3m(CLO_PATH)
opn_rec  = load_a3m(OPN_PATH)

print(f"Original MSA: {len(orig_rec)} sequences")
print(f"Closed500   : {len(clo_rec)} sequences")
print(f"Open100     : {len(opn_rec)} sequences")

# strip inserts (any lowercase letters are inserts vs the query)
orig_seqs = [strip_a3m_inserts(s) for _, s in orig_rec]
clo_seqs  = [strip_a3m_inserts(s) for _, s in clo_rec]
opn_seqs  = [strip_a3m_inserts(s) for _, s in opn_rec]

# distribution of stripped lengths
print("\nStripped sequence lengths:")
print(f"  orig: unique={set(len(s) for s in orig_seqs)}")
print(f"  clo : unique={set(len(s) for s in clo_seqs)}")
print(f"  opn : unique={set(len(s) for s in opn_seqs)}")

# queries: first record
print("\nQuery comparison (first sequence in each file):")
print(f"  orig[0] header : {orig_rec[0][0][:60]}")
print(f"  orig[0] seq    : {orig_seqs[0][:60]}...")
print(f"  clo[0]  header : {clo_rec[0][0][:60]}")
print(f"  clo[0]  seq    : {clo_seqs[0][:60]}...")
print(f"  opn[0]  header : {opn_rec[0][0][:60]}")
print(f"  opn[0]  seq    : {opn_seqs[0][:60]}...")
print(f"\norig[0] == clo[0]? {orig_seqs[0]==clo_seqs[0]}")
print(f"orig[0] == opn[0]? {orig_seqs[0]==opn_seqs[0]}")
print(f"clo[0]  == opn[0]? {clo_seqs[0]==opn_seqs[0]}")

# --- column-wise comparison at candidate columns ---
CANDIDATES = [136, 202, 172, 155, 192, 109, 135]    # 1-based residue numbers (closed-pushers + open-pusher 135)

print("\n" + "="*70)
print("Per-candidate column comparison (1-based residue indexing)")
print("="*70)
for res in CANDIDATES:
    col = res - 1   # 0-based
    if any(len(s)<=col for s in orig_seqs+clo_seqs+opn_seqs):
        print(f"\n  res {res}: SKIPPED (out of range in one of the files)")
        continue
    c_orig = col_dist(orig_seqs, col)
    c_clo  = col_dist(clo_seqs,  col)
    c_opn  = col_dist(opn_seqs,  col)
    f_orig = freqs(c_orig); f_clo = freqs(c_clo); f_opn = freqs(c_opn)
    top3 = lambda d: sorted(d.items(), key=lambda x:-x[1])[:3]
    print(f"\n  res {res}  query={orig_seqs[0][col]}")
    print(f"    orig    n={sum(c_orig.values()):4d}  top3: " +
          "  ".join(f"{a}:{p:.2f}" for a,p in top3(f_orig)))
    print(f"    Closed  n={sum(c_clo.values()):4d}  top3: " +
          "  ".join(f"{a}:{p:.2f}" for a,p in top3(f_clo)))
    print(f"    Open    n={sum(c_opn.values()):4d}  top3: " +
          "  ".join(f"{a}:{p:.2f}" for a,p in top3(f_opn)))
    js = js_divergence(f_clo, f_opn)
    print(f"    JS(Closed||Open) = {js:.3f}    "
          f"mode_closed={max(f_clo,key=f_clo.get) if f_clo else '?'}    "
          f"mode_open  ={max(f_opn,key=f_opn.get) if f_opn else '?'}    "
          f"discriminative={'YES' if max(f_clo,key=f_clo.get)!=max(f_opn,key=f_opn.get) else 'no'}")

# --- as a sanity-check, JS divergence at ALL columns, then rank ---
print("\n" + "="*70)
print("Global JS(Closed||Open) ranking (top 20 most discriminative columns)")
print("="*70)
L = 214
all_js = []
for col in range(L):
    fc = freqs(col_dist(clo_seqs, col))
    fo = freqs(col_dist(opn_seqs, col))
    if not fc or not fo: continue
    all_js.append((col+1, js_divergence(fc, fo), max(fc,key=fc.get), max(fo,key=fo.get)))
all_js.sort(key=lambda x:-x[1])
for res, js, mc, mo in all_js[:20]:
    print(f"  res {res:3d}  JS={js:.3f}  mode_closed={mc}  mode_open={mo}   {'switch!' if mc!=mo else '(same mode)'}")

# how many of our 7 candidates rank in top 20?
top20_res = {x[0] for x in all_js[:20]}
print(f"\nCandidates that land in JS top-20: {sorted(set(CANDIDATES) & top20_res)}")
print(f"Candidates' JS ranks: " + ", ".join(
    f"{res}=#{[r[0] for r in all_js].index(res)+1}({all_js[[r[0] for r in all_js].index(res)][1]:.2f})"
    for res in CANDIDATES if res in {r[0] for r in all_js}))
