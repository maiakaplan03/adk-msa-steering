#!/usr/bin/env python3
"""
02_delta_contacts.py — "2- Δcontactos".
Compares the contact sets of conformation A (1AKE, closed) and B (4AKE, open):
  lost   A−B : contacts only in 1AKE  -> delta_contacts_lost_A_minus_B.csv   (66, closed-specific)
  gained B−A : contacts only in 4AKE  -> delta_contacts_gained_B_minus_A.csv (57, open-specific)
  conserved  : in both                -> delta_contacts_conserved.csv        (423)
plus per-residue degree changes and the top gainers / losers with their partners.
Usage: python pipeline/02_delta_contacts.py --csv-a results/contacts/1AKE_residue_contacts_shadow.csv \
           --csv-b results/contacts/4AKE_residue_contacts_shadow.csv --out-dir results/contacts
"""
import argparse, os
import pandas as pd
from collections import Counter, defaultdict

# === Configurá acá los CSV salidos de tu script base ===
CSV_A = "1AKE_residue_contacts_shadow.csv"  # Conformación A: 1AKE (cerrada, con Ap5A)
CSV_B = "4AKE_residue_contacts_shadow.csv"  # Conformación B: 4AKE (abierta)
_ap = argparse.ArgumentParser()
_ap.add_argument("--csv-a", default=CSV_A)
_ap.add_argument("--csv-b", default=CSV_B)
_ap.add_argument("--out-dir", default=".")
_args = _ap.parse_args()
CSV_A, CSV_B = os.path.abspath(_args.csv_a), os.path.abspath(_args.csv_b)
os.makedirs(_args.out_dir, exist_ok=True)
os.chdir(_args.out_dir)   # all outputs are written next to each other, as in the notebook

# (Ajustá estos nombres si tus CSV usan otros)
COLS = dict(
    chain1="chain1", resn1="resn1", resi1="resi1",
    chain2="chain2", resn2="resn2", resi2="resi2",
    dist="dist_A"   # nombre de la columna de distancia en tus CSV
)

TOP_N = 10
TOP_PARTNER_LIST_LIMIT = 12  # para la impresión; en CSV se guarda completo

# -------------------- utilidades --------------------
def _pair_key(row):
    """((chain,resi),(chain,resi)) ordenado para que el contacto sea no dirigido."""
    a = (str(row[COLS["chain1"]]), int(row[COLS["resi1"]]))
    b = (str(row[COLS["chain2"]]), int(row[COLS["resi2"]]))
    return tuple(sorted([a, b]))

def _dedup_pairs(df):
    """Quita duplicados por par normalizado; conserva el de menor distancia si hay varios."""
    df = df.copy()
    df["_pair"] = df.apply(_pair_key, axis=1)
    if COLS["dist"] in df.columns:
        df = df.sort_values(COLS["dist"]).drop_duplicates(subset=["_pair"], keep="first")
    else:
        df = df.drop_duplicates(subset=["_pair"])
    return df.reset_index(drop=True)

def _lookup_info(df):
    """pair -> info (chains, resi, resn, dist) tomado de esa conformación."""
    look = {}
    for _, r in df.iterrows():
        p = r["_pair"]
        (c1, i1), (c2, i2) = p
        look[p] = dict(
            chain1=c1, resi1=int(i1), resn1=r.get(COLS["resn1"], pd.NA),
            chain2=c2, resi2=int(i2), resn2=r.get(COLS["resn2"], pd.NA),
            dist=r.get(COLS["dist"], pd.NA)
        )
    return look

def _pairs_to_df(pairs, info=None, dist_label="dist"):
    """Convierte set de pares a DataFrame legible, adjuntando info si está disponible."""
    rows = []
    for (a, b) in sorted(pairs):
        (c1, i1), (c2, i2) = a, b
        ii = {}
        if info and (a, b) in info: ii = info[(a, b)]
        elif info and (b, a) in info: ii = info[(b, a)]
        rows.append(dict(
            chain1=c1, resi1=int(i1), resn1=ii.get("resn1", pd.NA),
            chain2=c2, resi2=int(i2), resn2=ii.get("resn2", pd.NA),
            **{dist_label: ii.get("dist", pd.NA)}
        ))
    return pd.DataFrame(rows)

def degree_from_set(pairs):
    """Grado por residuo (número de contactos) a partir de un set de pares no dirigidos."""
    deg = Counter()
    for (a, b) in pairs:
        deg[a] += 1
        deg[b] += 1
    return deg

def adjacency_from_set(pairs):
    """Mapa de adyacencia: residuo -> set de socios conectados."""
    adj = defaultdict(set)
    for (a, b) in pairs:
        adj[a].add(b)
        adj[b].add(a)
    return adj

def build_resn_lookup(dfA, dfB):
    """(chain,resi) -> resn, recolectado de cualquiera de las tablas."""
    look = {}
    for df in (dfA, dfB):
        for _, r in df.iterrows():
            a = (str(r[COLS["chain1"]]), int(r[COLS["resi1"]]))
            b = (str(r[COLS["chain2"]]), int(r[COLS["resi2"]]))
            ra, rb = r.get(COLS["resn1"]), r.get(COLS["resn2"])
            if pd.notna(ra): look[a] = ra
            if pd.notna(rb): look[b] = rb
    return look

def fmt_res(res, resn_lookup):
    c, i = res
    resn = resn_lookup.get(res, "?")
    return f"{c}:{i}({resn})"

def partners_str(res, adj, resn_lookup, limit=None):
    parts = sorted(adj.get(res, []))
    display = parts[:limit] if (limit and len(parts) > limit) else parts
    s = ", ".join(fmt_res(p, resn_lookup) for p in display)
    if limit and len(parts) > limit:
        s += f" … (+{len(parts)-limit} más)"
    return s

# -------------------- carga y Δ --------------------
dfA = _dedup_pairs(pd.read_csv(CSV_A))
dfB = _dedup_pairs(pd.read_csv(CSV_B))

setA = set(dfA["_pair"])
setB = set(dfB["_pair"])

gained_pairs    = setB - setA   # aparecen en B
lost_pairs      = setA - setB   # desaparecen en B
conserved_pairs = setA & setB   # en ambas

infoA = _lookup_info(dfA)
infoB = _lookup_info(dfB)

# DataFrames con distancias
df_gained = _pairs_to_df(gained_pairs, info=infoB, dist_label="dist_B") \
    .sort_values(["chain1","resi1","chain2","resi2"]).reset_index(drop=True)
df_lost   = _pairs_to_df(lost_pairs,   info=infoA, dist_label="dist_A") \
    .sort_values(["chain1","resi1","chain2","resi2"]).reset_index(drop=True)

df_cons_A = _pairs_to_df(conserved_pairs, info=infoA, dist_label="dist_A")
df_cons_B = _pairs_to_df(conserved_pairs, info=infoB, dist_label="dist_B")
df_conserved = df_cons_A.merge(
    df_cons_B[["chain1","resi1","chain2","resi2","dist_B"]],
    on=["chain1","resi1","chain2","resi2"], how="left"
)
if "dist_A" in df_conserved and "dist_B" in df_conserved:
    df_conserved["delta_dist_B_minus_A"] = df_conserved["dist_B"] - df_conserved["dist_A"]
df_conserved = df_conserved.sort_values(["chain1","resi1","chain2","resi2"]).reset_index(drop=True)

# Guardar tablas base
df_gained.to_csv("delta_contacts_gained_B_minus_A.csv", index=False)
df_lost.to_csv("delta_contacts_lost_A_minus_B.csv", index=False)
df_conserved.to_csv("delta_contacts_conserved.csv", index=False)

# -------------------- degree y socios --------------------
degA = degree_from_set(setA)
degB = degree_from_set(setB)
adj_gained = adjacency_from_set(gained_pairs)  # socios nuevos por residuo
adj_lost   = adjacency_from_set(lost_pairs)    # socios perdidos por residuo
resn_lookup = build_resn_lookup(dfA, dfB)

# universo de residuos
all_res = sorted(set(list(degA.keys()) + list(degB.keys()) +
                     list(adj_gained.keys()) + list(adj_lost.keys())))

rows_deg = []
for res in all_res:
    c, i = res
    rows_deg.append({
        "chain": c,
        "resi": int(i),
        "resn": resn_lookup.get(res, pd.NA),
        "degree_A": degA.get(res, 0),
        "degree_B": degB.get(res, 0),
        "delta_degree_B_minus_A": degB.get(res, 0) - degA.get(res, 0),
        "gained_count": len(adj_gained.get(res, [])),
        "lost_count": len(adj_lost.get(res, [])),
        "gained_partners": ", ".join(fmt_res(p, resn_lookup) for p in sorted(adj_gained.get(res, []))),
        "lost_partners": ", ".join(fmt_res(p, resn_lookup) for p in sorted(adj_lost.get(res, []))),
    })
df_deg = pd.DataFrame(rows_deg).sort_values(["chain","resi"]).reset_index(drop=True)
df_deg.to_csv("delta_contacts_per_residue.csv", index=False)

# -------------------- RESUMEN RÁPIDO --------------------
print("\n=== RESUMEN ΔCONTACTOS ===")
print(f"Contactos A: {len(setA)} | Contactos B: {len(setB)}")
print(f"Ganados (B−A): {len(gained_pairs)} | Perdidos (A−B): {len(lost_pairs)} | Conservados: {len(conserved_pairs)}")

def _safe_stats(series):
    s = pd.to_numeric(series, errors="coerce").dropna()
    if len(s)==0: return "n/a"
    return f"n={len(s)}, mean={s.mean():.2f} Å, median={s.median():.2f}, min={s.min():.2f}, max={s.max():.2f}"

print("\n— Distancias en ganados (dist_B):", _safe_stats(df_gained.get("dist_B", pd.Series(dtype=float))))
print("— Distancias en perdidos (dist_A):", _safe_stats(df_lost.get("dist_A", pd.Series(dtype=float))))
if "delta_dist_B_minus_A" in df_conserved.columns:
    print("— Δdist conservados (B−A):", _safe_stats(df_conserved["delta_dist_B_minus_A"]))

# -------------------- TOPS con socios --------------------
top_gain = df_deg.sort_values("delta_degree_B_minus_A", ascending=False).head(TOP_N).copy()
top_loss = df_deg.sort_values("delta_degree_B_minus_A", ascending=True).head(TOP_N).copy()

# versiones “cortas” para impresión
top_gain["gained_partners_preview"] = [
    partners_str((r["chain"], r["resi"]), adj_gained, resn_lookup, limit=TOP_PARTNER_LIST_LIMIT)
    for _, r in top_gain.iterrows()
]
top_loss["lost_partners_preview"] = [
    partners_str((r["chain"], r["resi"]), adj_lost, resn_lookup, limit=TOP_PARTNER_LIST_LIMIT)
    for _, r in top_loss.iterrows()
]

print(f"\n— Top {TOP_N} residuos que MÁS ganan contactos (Δgrado):")
print(top_gain[["chain","resi","resn","degree_A","degree_B","delta_degree_B_minus_A",
                "gained_count","lost_count","gained_partners_preview"]].to_string(index=False))

print(f"\n— Top {TOP_N} residuos que MÁS pierden contactos (Δgrado):")
print(top_loss[["chain","resi","resn","degree_A","degree_B","delta_degree_B_minus_A",
                "gained_count","lost_count","lost_partners_preview"]].to_string(index=False))

# Guardar los tops completos con listas completas de socios
top_gain.to_csv("delta_top_gainers_with_partners.csv", index=False)
top_loss.to_csv("delta_top_losers_with_partners.csv", index=False)

print("\nArchivos escritos:")
print(" - delta_contacts_gained_B_minus_A.csv")
print(" - delta_contacts_lost_A_minus_B.csv")
print(" - delta_contacts_conserved.csv")
print(" - delta_contacts_per_residue.csv  (incluye resn + socios ganados/perdidos)")
print(" - delta_top_gainers_with_partners.csv")
print(" - delta_top_losers_with_partners.csv")
