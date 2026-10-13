#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05a_qa_qb.py — "5.1 análisis de estructuras" (QA vs QB map).

Usage:
  python pipeline/05a_qa_qb.py --structures-dir all_structures --results-dir results_contactos_shadow
The structures folder must also contain the two references, 1AKE_A.pdb and 4AKE_A.pdb.

Pipeline automático para calcular contactos "Shadow" (átomos pesados + oclusión)
a partir de múltiples estructuras .cif/.pdb, comparar contactos entre
dos conformaciones de referencia (A y B), y generar un gráfico QA vs QB.

Criterio Shadow:
- CUTOFF C (Å) para vecinos átomo–átomo.
- Se descartan si hay un átomo intermedio cuya distancia al segmento <= S (Å).
- Se excluyen vecinos 1–2 y 1–3 en la misma cadena (MIN_SEQ_SEP=3).
- Se trabaja a nivel residuo–residuo (pares únicos normalizados).
"""

# ===========================================================
# ✅ CONFIGURACIÓN
# ===========================================================
import argparse
import matplotlib
matplotlib.use("Agg")
from Bio.PDB import MMCIFParser, PDBParser, is_aa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.spatial import cKDTree as KDTree

# --- Rutas / carpetas ---
_ap = argparse.ArgumentParser()
_ap.add_argument("--structures-dir", default="all_structures")
_ap.add_argument("--results-dir", default="results_contactos_shadow")
_ap.add_argument("--reference-a", default="1AKE_A.pdb")
_ap.add_argument("--reference-b", default="4AKE_A.pdb")
_args = _ap.parse_args()
STRUCTURES_DIR = _args.structures_dir
RESULTS_DIR = _args.results_dir
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- Modelos de referencia (cerrada / abierta) ---
REFERENCE_A = _args.reference_a   # 1AKE_A.pdb (cerrada)
REFERENCE_B = _args.reference_b   # 4AKE_A.pdb (abierta)

# --- Parámetros Shadow ---
CUTOFF = 5.2          # C (Å) - vecindad átomo–átomo
SHADOW_S = 1.3        # S (Å) - radio de "sombreado" para oclusión
MIN_SEQ_SEP = 3       # excluir 1–2 y 1–3 (|Δseq| <= 2)
CHAIN_ID = None       # "A" para restringir, None = todas
ALLOW_INTERCHAINS = False
EXCLUDE_WATER = True

# ===========================================================
# ✅ FUNCIONES DE E/S
# ===========================================================
def load_structure(path):
    """Carga estructura con Bio.Python, soporta .cif y .pdb. Devuelve el primer modelo."""
    ext = os.path.splitext(path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("model", path)
        return structure[0]
    except Exception as e:
        print(f"Error loading structure {path}: {e}")
        return None

def is_heavy(atom):
    """Excluye H/D/T si aparecieran."""
    el = (atom.element or "").upper()
    if not el:
        # fallback por nombre
        name = atom.get_name().strip()
        el = "".join([c for c in name if c.isalpha()])[:1].upper()
    return el not in {"H", "D", "T"}

# ===========================================================
# ✅ UTILIDADES GEOMÉTRICAS (Shadow)
# ===========================================================
def point_segment_distance(p, a, b):
    """Distancia de p al segmento ab y parámetro t (proyección)."""
    ab = b - a
    ap = p - a
    ab2 = float(np.dot(ab, ab))
    if ab2 == 0.0:
        return float(np.linalg.norm(ap)), 0.0
    t = max(0.0, min(1.0, float(np.dot(ap, ab) / ab2)))
    proj = a + t * ab
    return float(np.linalg.norm(p - proj)), t

# ===========================================================
# ✅ EXTRACCIÓN DE ÁTOMOS Y ORDEN SECUENCIAL
# ===========================================================
def chain_ok(chain_id):
    return True if CHAIN_ID is None else (chain_id == CHAIN_ID)

def harvest_heavy_atoms(model):
    """
    Recolecta átomos pesados por residuo e índice secuencial por cadena.
    Devuelve:
    - atoms: lista de dicts con info y coords
    - res_order: (chain, resseq, icode) -> índice secuencial
    """
    # construir orden de residuos por cadena
    res_order = {}
    for ch in model:
        if not chain_ok(ch.id):
            continue
        residues = [r for r in ch if (not EXCLUDE_WATER or r.id[0] != "W") and is_aa(r, standard=True)]
        residues.sort(key=lambda r: (r.id[1], r.id[2] or " "))
        for i, r in enumerate(residues):
            key = (ch.id, r.id[1], r.id[2] or " ")
            res_order[key] = i

    atoms = []
    idx = 0
    for ch in model:
        if not chain_ok(ch.id):
            continue
        for r in ch:
            if EXCLUDE_WATER and r.id[0] == "W":
                continue
            if not is_aa(r, standard=True):
                continue
            key = (ch.id, r.id[1], r.id[2] or " ")
            seq_index = res_order.get(key, None)
            for a in r:
                if not is_heavy(a):
                    continue
                atoms.append(dict(
                    idx=idx,
                    chain=ch.id,
                    resname=r.get_resname(),
                    resseq=r.id[1],
                    icode=r.id[2] or " ",
                    atom=a.get_name().strip(),
                    is_backbone=a.get_name().strip() in {"N","CA","C","O","OXT"},
                    coord=a.get_coord().astype(float),
                    seq_index=seq_index
                ))
                idx += 1
    return atoms, res_order

def normalize_res_pair(p):
    """Normaliza un par de residuos a una tupla ordenada estable para usar en sets/dicts."""
    (ci, ri, ii), (cj, rj, ij) = p
    a = (ci, int(ri), str(ii))
    b = (cj, int(rj), str(ij))
    return tuple(sorted([a, b]))

# ===========================================================
# ✅ CONTACTOS SHADOW (RESIDUO–RESIDUO)
# ===========================================================
def shadow_contacts_residue(model, C=CUTOFF, S=SHADOW_S,
                            min_seq_sep=MIN_SEQ_SEP,
                            allow_interchains=ALLOW_INTERCHAINS):
    """
    Calcula contactos residuo–residuo con criterio Shadow.
    Devuelve: set de pares ((chain_i, res_i, icode_i), (chain_j, res_j, icode_j))
    """
    atoms, _ = harvest_heavy_atoms(model)
    if not atoms:
        return set()

    coords = np.vstack([a["coord"] for a in atoms])
    tree = KDTree(coords)

    res_pairs = set()

    for i, ai in enumerate(atoms):
        # buscar vecinos candidatos en radio C
        neigh = tree.query_ball_point(ai["coord"], C)
        for j in neigh:
            if j <= i:
                continue
            aj = atoms[j]

            # filtro inter-cadenas
            if not allow_interchains and ai["chain"] != aj["chain"]:
                continue

            # excluir 1–2 y 1–3 (misma cadena)
            if ai["chain"] == aj["chain"] and ai["seq_index"] is not None and aj["seq_index"] is not None:
                if abs(ai["seq_index"] - aj["seq_index"]) <= (min_seq_sep - 1):
                    continue

            # distancia real
            d = float(np.linalg.norm(ai["coord"] - aj["coord"]))
            if d > C or d <= 1e-6:
                continue

            # candidatos para oclusión: esfera alrededor del punto medio
            midpoint = 0.5 * (ai["coord"] + aj["coord"])
            search_r = 0.5 * d + S
            cand = tree.query_ball_point(midpoint, search_r)

            # Ajuste conservador: si son casi secuenciales (p. ej., 1–4), reducir S
            S_eff = S
            if ai["chain"] == aj["chain"] and ai["seq_index"] is not None and aj["seq_index"] is not None:
                if abs(ai["seq_index"] - aj["seq_index"]) == min_seq_sep:
                    S_eff = min(S, 0.5)

            # test de sombra: existe k que ocluye?
            A = ai["coord"]; B = aj["coord"]
            shadowed = False
            for k in cand:
                if k == i or k == j:
                    continue
                P = coords[int(k)]
                # prune rápido por cercanía a extremos
                if np.linalg.norm(P - A) >= d and np.linalg.norm(P - B) >= d:
                    continue
                dist, t = point_segment_distance(P, A, B)
                if dist <= S_eff and 0.0 <= t <= 1.0:
                    shadowed = True
                    break
            if shadowed:
                continue

            # si pasó el filtro, agregar par residuo–residuo
            pair = normalize_res_pair((
                (ai["chain"], ai["resseq"], ai["icode"]),
                (aj["chain"], aj["resseq"], aj["icode"])
            ))
            res_pairs.add(pair)

    return res_pairs

def compute_contact_sets(structure_path):
    model = load_structure(structure_path)
    if model is None:
        return set()
    return shadow_contacts_residue(model, C=CUTOFF, S=SHADOW_S,
                                   min_seq_sep=MIN_SEQ_SEP,
                                   allow_interchains=ALLOW_INTERCHAINS)

# ===========================================================
# ✅ CÁLCULO DE CONTACTOS EN LOTE
# ===========================================================
all_files = [f for f in os.listdir(STRUCTURES_DIR) if f.lower().endswith((".cif", ".pdb"))]
print(f"\nSe encontraron {len(all_files)} estructuras:")
for f in all_files:
    print(" -", f)

contacts_dict = {}
for i, fname in enumerate(sorted(all_files), 1):
    path = os.path.join(STRUCTURES_DIR, fname)
    print(f"\n[{i}/{len(all_files)}] Procesando {fname} ...", end=" ")
    contacts = compute_contact_sets(path)
    print(f"{len(contacts)} contactos residuo–residuo (Shadow).")
    contacts_dict[fname] = contacts

# ===========================================================
# ✅ DELTA-CONTACTOS Y QA/QB
# ===========================================================
if REFERENCE_A not in contacts_dict or REFERENCE_B not in contacts_dict:
    raise KeyError(
        f"Las referencias no están entre las estructuras cargadas.\n"
        f"REFERENCE_A={REFERENCE_A}\nREFERENCE_B={REFERENCE_B}\n"
        f"Claves disponibles: {list(contacts_dict.keys())}"
    )

A = contacts_dict[REFERENCE_A]
B = contacts_dict[REFERENCE_B]
A_only = A - B
B_only = B - A
print(f"\nContactos exclusivos de A: {len(A_only)} | exclusivos de B: {len(B_only)}")

# Guardar conjuntos A_only/B_only para inspección
def pairs_to_df(pairs):
    rows = []
    for (ci, ri, ii), (cj, rj, ij) in pairs:
        rows.append(dict(chain_i=ci, res_i=ri, icode_i=ii, chain_j=cj, res_j=rj, icode_j=ij))
    return pd.DataFrame(rows)

pairs_to_df(A_only).to_csv(os.path.join(RESULTS_DIR, "A_only_pairs.csv"), index=False)
pairs_to_df(B_only).to_csv(os.path.join(RESULTS_DIR, "B_only_pairs.csv"), index=False)

rows = []
for fname, contacts in contacts_dict.items():
    QA = len(contacts & A_only)
    QB = len(contacts & B_only)
    QA_frac = QA / len(A_only) if len(A_only) else float("nan")
    QB_frac = QB / len(B_only) if len(B_only) else float("nan")
    rows.append(dict(
        model_label=os.path.splitext(fname)[0],
        qa=QA, qb=QB,
        qa_frac=QA_frac, qb_frac=QB_frac,
        total_contacts=len(contacts)
    ))

summary = pd.DataFrame(rows).sort_values(["qb", "qa"], ascending=[False, True])
summary_path = os.path.join(RESULTS_DIR, "QA_QB_summary_shadow.csv")
summary.to_csv(summary_path, index=False)

print("\nResumen QA/QB por modelo:")
print(summary.to_string(index=False))

# ===========================================================
# ✅ GRÁFICO QA vs QB (solo 1AKE_A y 4AKE_A con nombre; resto numeradas)
# ===========================================================
import os

plt.figure(figsize=(6.8, 6.8))
plt.scatter(summary["qa"], summary["qb"])

# Identificadores base (sin extensión) de las referencias
ref_a = os.path.splitext(REFERENCE_A)[0]
ref_b = os.path.splitext(REFERENCE_B)[0]

# Asignar números SOLO a los no-referencia en el orden en que aparecen en el scatter (summary)
num_map = {}   # model_label -> número
counter = 1
for _, r in summary.iterrows():
    label = r["model_label"]
    if label not in {ref_a, ref_b}:
        num_map[label] = counter
        counter += 1

# Etiquetado: refs por nombre, resto por número
for _, r in summary.iterrows():
    x, y = float(r["qa"]), float(r["qb"])
    label = r["model_label"]
    if label in {ref_a, ref_b}:
        plt.text(x, y, label, fontsize=9)
    else:
        plt.text(x, y, str(num_map[label]), fontsize=8)

plt.xlabel("QA: contacts exclusive to A (1AKE) present in the model (Shadow)")
plt.ylabel("QB: contacts exclusive to B (4AKE) present in the model (Shadow)")
plt.title("QA vs QB map (Shadow criterion)")

# Puntos “ideal” A y B
plt.scatter([len(A_only)], [0], marker="x")
plt.text(len(A_only), 0, "A (ideal)", fontsize=8, ha="left", va="bottom")
plt.scatter([0], [len(B_only)], marker="x")
plt.text(0, len(B_only), "B (ideal)", fontsize=8, ha="left", va="bottom")

plt.tight_layout()
plot_path = os.path.join(RESULTS_DIR, "QA_vs_QB_shadow.png")
plt.savefig(plot_path, dpi=150)
plt.show()

# Guardar leyenda de números → modelos (excluye referencias)
legend_rows = [
    {"number": n, "model_label": lbl}
    for lbl, n in sorted(num_map.items(), key=lambda x: x[1])
]
legend_df = pd.DataFrame(legend_rows)
legend_path = os.path.join(RESULTS_DIR, "legend_QA_QB_numbers.csv")
legend_df.to_csv(legend_path, index=False)

print("Etiquetado OK. Archivos guardados:")
print(" -", plot_path)
print(" -", legend_path, "(número -> modelo para los no-referencia)")
