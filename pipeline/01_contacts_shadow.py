#!/usr/bin/env python3
"""
01_contacts_shadow.py — "1- cálculo de contactos".
Residue-residue contacts of one structure with the Shadow criterion:
heavy atoms within 5.2 Å (C), discarded if another atom lies within 1.3 Å (S) of the
segment joining them, and |Δseq| >= 3.
Usage:
  python pipeline/01_contacts_shadow.py --pdb pdb_chainA_protein_only/1AKE_chainA_protein_only.pdb \
         --out-contacts 1AKE_residue_contacts_shadow.csv --out-matrix contact_matrix_1AKE_chainA_shadow.csv
Outputs for 1AKE and 4AKE are in results/contacts/.
"""
import argparse
from Bio.PDB import PDBParser, is_aa
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree as KDTree

# =========================
# CONFIG
# =========================
PDB_PATH = "/content/pdb_chainA_protein_only/1AKE_chainA_protein_only.pdb"   # cambia esto

CUTOFF = 5.2          # C (Å): vecindad átomo–átomo
SHADOW_S = 1.3        # S (Å): radio de oclusión alrededor del segmento
MIN_SEQ_SEP = 3       # excluye 1–2 y 1–3 (|Δseq| <= 2) en la misma cadena

TARGET_CHAIN = "A"    # para la matriz (puede ser None si no querés matriz)
CHAIN_ID = None       # None = usar todas las cadenas; "A" = restringir análisis a una cadena
ALLOW_INTERCHAINS = False  # True = permite contactos entre cadenas (útil para cristal)
EXCLUDE_WATER = True       # True = excluye HOH/WAT
STANDARD_AA_ONLY = True    # True = solo aa estándar

OUT_CONTACTS = "1AKE_residue_contacts_shadow.csv"
OUT_MATRIX = "contact_matrix_1AKE_chainA_shadow.csv"

# =========================
# HELPERS
# =========================
def chain_ok(cid: str) -> bool:
    return True if CHAIN_ID is None else (cid == CHAIN_ID)

def is_heavy(atom) -> bool:
    """Excluye H/D/T si aparecieran; fallback por nombre si atom.element no está."""
    el = (atom.element or "").upper()
    if not el:
        name = atom.get_name().strip()
        el = "".join([c for c in name if c.isalpha()])[:1].upper()
    return el not in {"H", "D", "T"}

def point_segment_distance(p, a, b):
    """Distancia de p al segmento ab y parámetro t de proyección."""
    ab = b - a
    ap = p - a
    ab2 = float(np.dot(ab, ab))
    if ab2 == 0.0:
        return float(np.linalg.norm(ap)), 0.0
    t = float(np.dot(ap, ab) / ab2)
    t = max(0.0, min(1.0, t))
    proj = a + t * ab
    return float(np.linalg.norm(p - proj)), t

def resid_key(chain_id, res):
    """Clave estable de residuo: (chain, resseq, icode)."""
    return (chain_id, int(res.id[1]), str(res.id[2] or " "))

def normalize_pair(a, b):
    """Orden estable para pares de residuos."""
    return tuple(sorted([a, b]))

def harvest_heavy_atoms(model):
    """
    Devuelve:
      atoms: lista de dicts con coords y metadata
      res_order: dict (chain, resseq, icode)-> índice secuencial por cadena (para MIN_SEQ_SEP)
      res_info: dict (chain, resseq, icode)-> resname
    """
    res_order = {}
    res_info = {}

    # Orden secuencial por cadena
    for ch in model:
        if not chain_ok(ch.id):
            continue
        residues = []
        for r in ch:
            if EXCLUDE_WATER and r.id[0] in {"W"}:
                continue
            if STANDARD_AA_ONLY:
                if not is_aa(r, standard=True):
                    continue
            else:
                if not is_aa(r):
                    continue
            residues.append(r)

        residues.sort(key=lambda r: (r.id[1], r.id[2] or " "))
        for i, r in enumerate(residues):
            k = resid_key(ch.id, r)
            res_order[k] = i
            res_info[k] = r.get_resname()

    # Átomos
    atoms = []
    idx = 0
    for ch in model:
        if not chain_ok(ch.id):
            continue
        for r in ch:
            if EXCLUDE_WATER and r.id[0] in {"W"}:
                continue
            if STANDARD_AA_ONLY:
                if not is_aa(r, standard=True):
                    continue
            else:
                if not is_aa(r):
                    continue

            rk = resid_key(ch.id, r)
            seq_index = res_order.get(rk, None)

            for a in r:
                if not is_heavy(a):
                    continue
                atoms.append(
                    dict(
                        idx=idx,
                        chain=ch.id,
                        reskey=rk,
                        resname=r.get_resname(),
                        atom=a.get_name().strip(),
                        coord=a.get_coord().astype(float),
                        seq_index=seq_index,
                    )
                )
                idx += 1

    return atoms, res_order, res_info

# =========================
# SHADOW CONTACTS (single structure)
# =========================
def shadow_contacts_residue_single(model):
    atoms, _, res_info = harvest_heavy_atoms(model)
    if not atoms:
        return {}, res_info

    coords = np.vstack([a["coord"] for a in atoms])
    tree = KDTree(coords)

    # Par residuo-residuo -> distancia mínima (entre átomos pesados NO-sombreados)
    pair_min_d = {}

    for i, ai in enumerate(atoms):
        neigh = tree.query_ball_point(ai["coord"], CUTOFF)
        for j in neigh:
            if j <= i:
                continue
            aj = atoms[j]

            # Filtro inter-cadenas
            if (not ALLOW_INTERCHAINS) and (ai["chain"] != aj["chain"]):
                continue

            # Excluir 1–2 y 1–3 (misma cadena)
            if ai["chain"] == aj["chain"] and ai["seq_index"] is not None and aj["seq_index"] is not None:
                if abs(ai["seq_index"] - aj["seq_index"]) <= (MIN_SEQ_SEP - 1):
                    continue

            A = ai["coord"]
            B = aj["coord"]
            d = float(np.linalg.norm(A - B))
            if d > CUTOFF or d <= 1e-6:
                continue

            # Candidatos para oclusión alrededor del punto medio
            midpoint = 0.5 * (A + B)
            search_r = 0.5 * d + SHADOW_S
            cand = tree.query_ball_point(midpoint, search_r)

            # Ajuste conservador como tu pipeline (si Δseq == MIN_SEQ_SEP)
            S_eff = SHADOW_S
            if ai["chain"] == aj["chain"] and ai["seq_index"] is not None and aj["seq_index"] is not None:
                if abs(ai["seq_index"] - aj["seq_index"]) == MIN_SEQ_SEP:
                    S_eff = min(SHADOW_S, 0.5)

            shadowed = False
            for k in cand:
                if k == i or k == j:
                    continue
                P = coords[int(k)]
                # prune rápido: si está lejos de ambos extremos
                if np.linalg.norm(P - A) >= d and np.linalg.norm(P - B) >= d:
                    continue
                dist, t = point_segment_distance(P, A, B)
                if dist <= S_eff and 0.0 <= t <= 1.0:
                    shadowed = True
                    break
            if shadowed:
                continue

            # Registrar contacto a nivel residuo-residuo
            rpi = ai["reskey"]
            rpj = aj["reskey"]
            pair = normalize_pair(rpi, rpj)

            prev = pair_min_d.get(pair, None)
            if prev is None or d < prev:
                pair_min_d[pair] = d

    return pair_min_d, res_info

def build_contact_df(pair_min_d, res_info):
    rows = []
    for (ci, ri, ii), (cj, rj, ij), in pair_min_d.keys():
        d = pair_min_d[((ci, ri, ii), (cj, rj, ij))]
        rows.append(
            dict(
                chain1=ci,
                resn1=res_info.get((ci, ri, ii), ""),
                resi1=int(ri),
                icode1=str(ii),
                chain2=cj,
                resn2=res_info.get((cj, rj, ij), ""),
                resi2=int(rj),
                icode2=str(ij),
                dist_A=float(d),  # “distancia mínima” entre átomos pesados NO-sombreados
            )
        )
    df = pd.DataFrame(rows).sort_values(["chain1", "resi1", "icode1", "chain2", "resi2", "icode2"])
    return df

def build_chain_matrix(df_contacts, chain_id, res_info):
    # residuos de esa cadena (únicos) ordenados por (resseq, icode)
    keys = sorted(
        [k for k in res_info.keys() if k[0] == chain_id],
        key=lambda x: (x[1], x[2])
    )
    key_to_idx = {k: i for i, k in enumerate(keys)}
    n = len(keys)
    M = np.zeros((n, n), dtype=int)

    # Solo contactos intra-cadena en esa cadena
    for _, r in df_contacts.iterrows():
        if r["chain1"] == chain_id and r["chain2"] == chain_id:
            k1 = (chain_id, int(r["resi1"]), str(r["icode1"]))
            k2 = (chain_id, int(r["resi2"]), str(r["icode2"]))
            if k1 in key_to_idx and k2 in key_to_idx:
                i = key_to_idx[k1]
                j = key_to_idx[k2]
                M[i, j] = 1
                M[j, i] = 1

    # opcional: devolver también etiquetas
    labels = [f"{k[1]}{'' if k[2].strip()=='' else k[2].strip()}" for k in keys]
    return M, labels

# =========================
# RUN
# =========================
_ap = argparse.ArgumentParser()
_ap.add_argument("--pdb", default=PDB_PATH)
_ap.add_argument("--out-contacts", default=OUT_CONTACTS)
_ap.add_argument("--out-matrix", default=OUT_MATRIX)
_args = _ap.parse_args()
PDB_PATH, OUT_CONTACTS, OUT_MATRIX = _args.pdb, _args.out_contacts, _args.out_matrix

parser = PDBParser(QUIET=True)
structure = parser.get_structure("prot", PDB_PATH)
model = structure[0]

pair_min_d, res_info = shadow_contacts_residue_single(model)
df = build_contact_df(pair_min_d, res_info)
df.to_csv(OUT_CONTACTS, index=False)
print(f"OK: {len(df)} contactos (Shadow) -> {OUT_CONTACTS}")

if TARGET_CHAIN is not None:
    M, labels = build_chain_matrix(df, TARGET_CHAIN, res_info)
    # guardo matriz como CSV simple (como tu script), sin labels
    np.savetxt(OUT_MATRIX, M, fmt="%d", delimiter=",")
    print(f"OK: matriz de contactos (Shadow) cadena {TARGET_CHAIN} -> {OUT_MATRIX}")
    print(f"   tamaño: {M.shape[0]} x {M.shape[1]} (residuos en cadena {TARGET_CHAIN})")
