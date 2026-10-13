#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
05b_pca.py — "5.2 análisis de estructuras - PCA".

Usage:
  python pipeline/05b_pca.py --structures-dir all_structures --results-dir results_pca \
         [--noise-pc1-min 40 --noise-pc1-max 150]

Step 1: Kabsch-superpose every structure on 4AKE_A (Cα), PCA of the flattened Cα coordinates.
Step 2: mark the structures whose PC1 falls in [noise-pc1-min, noise-pc1-max] as "ruido"
        (the separate cluster in the first PCA) and redo the PCA without them.
        The window is specific to one run: look at the first plot before reusing it.

PCA estructural (Cα) con alineamiento rígido previo (Kabsch) usando Bio.PDB.Superimposer.
Lee .pdb/.cif desde STRUCTURES_DIR, alinea todo a la referencia y proyecta PC1–PC2.
Guarda además el RMSD de superposición a la referencia.
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, MMCIFParser, is_aa, Superimposer

# ===========================================================
# ✅ CONFIGURACIÓN
# ===========================================================
_ap = argparse.ArgumentParser()
_ap.add_argument("--structures-dir", default="all_structures")
_ap.add_argument("--results-dir", default="results_pca")
_ap.add_argument("--reference", default="4AKE_A.pdb")
_ap.add_argument("--noise-pc1-min", type=float, default=40.0)
_ap.add_argument("--noise-pc1-max", type=float, default=150.0)
_args = _ap.parse_args()
STRUCTURES_DIR = _args.structures_dir
RESULTS_DIR    = _args.results_dir
os.makedirs(RESULTS_DIR, exist_ok=True)

# Antes: CHAIN_ID = "A"
REFERENCE = _args.reference   # 4AKE_A.pdb por defecto; podés cambiar la referencia

# ===========================================================
# ✅ FUNCIONES AUXILIARES
# ===========================================================
def load_structure_generic(path):
    """Carga .pdb o .cif según extensión y retorna el primer modelo."""
    ext = os.path.splitext(path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure("model", path)
    return structure[0]

def get_main_chain(model):
    """
    Devuelve la cadena proteica 'principal'.
    Como todas tus estructuras tienen una sola cadena, en la práctica será esa,
    pero por robustez se elige la cadena con más residuos estándar.
    """
    chains = list(model.get_chains())
    if not chains:
        return None

    def count_std_res(chain):
        return sum(1 for r in chain if r.id[0] == " " and is_aa(r, standard=True))

    return max(chains, key=count_std_res)

def get_ca_atoms(model):
    """
    Lista de átomos CA (orden de recorrido) para la cadena principal.
    No depende del nombre de la cadena (A, 4AKE, etc.).
    """
    atoms = []
    chain = get_main_chain(model)
    if chain is None:
        return atoms
    for res in chain:
        if res.id[0] == " " and "CA" in res and is_aa(res, standard=True):
            atoms.append(res["CA"])
    return atoms

def coords_from_atoms(atoms):
    """np.array (N,3) a partir de lista de átomos."""
    return np.array([a.get_coord() for a in atoms], dtype=float)

def superimpose_to_ref(mobile_model, ref_model):
    """
    Alinea 'mobile_model' a 'ref_model' con Kabsch usando Cα.
    Modifica el 'mobile_model' in-place. Devuelve RMSD (float) o None si no pudo.
    """
    ref_atoms = get_ca_atoms(ref_model)
    mob_atoms = get_ca_atoms(mobile_model)
    if len(ref_atoms) == 0 or len(mob_atoms) == 0 or len(ref_atoms) != len(mob_atoms):
        return None
    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)
    sup.apply(mobile_model.get_atoms())
    return float(sup.rms)

# ===========================================================
# ✅ CARGA REFERENCIA
# ===========================================================
files = [f for f in os.listdir(STRUCTURES_DIR) if f.lower().endswith((".pdb", ".cif"))]
assert len(files) > 0, "No se encontraron estructuras en STRUCTURES_DIR"

ref_path = os.path.join(STRUCTURES_DIR, REFERENCE) if REFERENCE in files else os.path.join(STRUCTURES_DIR, files[0])
ref_model = load_structure_generic(ref_path)
ref_cas = get_ca_atoms(ref_model)
n_ref = len(ref_cas)
if n_ref == 0:
    raise ValueError(
        f"La referencia {os.path.basename(ref_path)} no tiene Cα estándar en su cadena proteica principal."
    )

print(f"Referencia: {os.path.basename(ref_path)} | Cα: {n_ref}")

# ===========================================================
# ✅ ALINEAMIENTO + EXTRACCIÓN DE COORDENADAS
# ===========================================================
X = []            # matriz (n_structures × 3N)
labels = []       # nombres para cada estructura
rmsds = []        # RMSD de alineamiento vs referencia

print("\nAlineando y extrayendo Cα:")
for f in files:
    path = os.path.join(STRUCTURES_DIR, f)
    try:
        model = load_structure_generic(path)

        # skip si no coincide el número de Cα
        if len(get_ca_atoms(model)) != n_ref:
            print(f" ⚠️ {f}: Nº de Cα distinto a la referencia. Se omite.")
            continue

        rms = superimpose_to_ref(model, ref_model)
        if rms is None:
            print(f" ⚠️ {f}: No se pudo alinear (longitudes distintas). Se omite.")
            continue

        ca_coords = coords_from_atoms(get_ca_atoms(model))  # ya alineadas
        X.append(ca_coords.flatten())
        labels.append(os.path.splitext(f)[0])
        rmsds.append(rms)
        print(f" - {f}: alineado (RMSD={rms:.3f} Å)")
    except Exception as e:
        print(f" ⚠️ {f}: error {e}. Se omite.")

if len(X) < 2:
    raise RuntimeError("Necesito al menos 2 estructuras alineadas para hacer PCA.")

X = np.vstack(X)  # (m × 3N)

# ===========================================================
# ✅ PCA (hasta 3 componentes)
# ===========================================================
n_comps = min(3, X.shape[0], X.shape[1])  # no más de 3, ni más de muestras o dims

try:
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_comps)
    proj = pca.fit_transform(X)  # (m × n_comps)
    var_ratio = pca.explained_variance_ratio_
except Exception:
    # Fallback SVD
    Xc = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    n_comps = min(n_comps, U.shape[1])
    proj = U[:, :n_comps] * S[:n_comps]
    var = (S**2) / (X.shape[0]-1)
    var_ratio = (var / var.sum())[:n_comps]

# Aseguramos columnas PC1..PC3 si existen
cols = {"model": labels, "RMSD_to_ref": rmsds}
for i in range(n_comps):
    cols[f"PC{i+1}"] = proj[:, i]

df_proj = pd.DataFrame(cols).sort_values("model").reset_index(drop=True)
df_proj.to_csv(os.path.join(RESULTS_DIR, "pca_projection_kabsch.csv"), index=False)

# === Etiquetado: 1AKE_A y 4AKE_A por nombre; el resto numerado ===
SPECIAL = {"1AKE_A", "4AKE_A"}   # sin extensión, coincide con df_proj["model"]

display_labels = []
number_map = []   # para guardar el mapeo número -> modelo
counter = 1
for name in df_proj["model"]:
    if name in SPECIAL:
        display_labels.append(name)
    else:
        display_labels.append(str(counter))
        number_map.append((counter, name))
        counter += 1

df_proj["display"] = display_labels

# Guardar mapeo números → nombres (solo los no especiales)
if number_map:
    map_df = pd.DataFrame(number_map, columns=["label_number", "model"])
    map_df.to_csv(os.path.join(RESULTS_DIR, "pca_label_map.csv"), index=False)

# Guardamos varianzas (si existen)
var1 = float(var_ratio[0]) if n_comps >= 1 else None
var2 = float(var_ratio[1]) if n_comps >= 2 else None
var3 = float(var_ratio[2]) if n_comps >= 3 else None

print("\nVarianza explicada:")
for i in range(n_comps):
    print(f"PC{i+1}: {var_ratio[i]*100:.2f}%")

# ===========================================================
# ✅ GRÁFICO PC1–PC2 (con nombres para 1AKE_A/4AKE_A y números para el resto)
# ===========================================================
plt.figure(figsize=(7, 7))

is_special = df_proj["model"].isin(SPECIAL)
# puntos no especiales
plt.scatter(df_proj.loc[~is_special, "PC1"], df_proj.loc[~is_special, "PC2"], alpha=0.8)
# puntos especiales resaltados
plt.scatter(df_proj.loc[is_special, "PC1"], df_proj.loc[is_special, "PC2"], s=80, marker="*", zorder=3)

# Anotaciones
for _, r in df_proj.iterrows():
    txt = r["display"]
    kws = dict(xytext=(3, 3), textcoords="offset points", fontsize=8)
    if r["model"] in SPECIAL:
        kws.update(fontweight="bold")
    plt.annotate(txt, (r["PC1"], r["PC2"]), **kws)

plt.xlabel(f"PC1 ({var1*100:.1f}% var)")
plt.ylabel(f"PC2 ({var2*100:.1f}% var)")
plt.title(f"Cα PCA with Kabsch superposition (ref: {os.path.basename(ref_path)})")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(RESULTS_DIR, "pca_projection_kabsch.png"), dpi=150)
plt.show()

print(" - pca_label_map.csv")


# ===========================================================
# 🔎 Filtrado por "ruido" en PC1 y PCA solo con las no-ruido
#    - Ruido:    noise_pc1_min <= PC1 <= noise_pc1_max (default 40..150)
#    - Exporta listas y rehace PCA con las no-ruido
#    - Grafica solo PC1 vs PC2
# ===========================================================

# --- 1) Marcar ruido según PC1 ---
mask_ruido = df_proj["PC1"].between(_args.noise_pc1_min, _args.noise_pc1_max, inclusive="both")
model_ruido    = df_proj.loc[mask_ruido, "model"].tolist()
model_no_ruido = df_proj.loc[~mask_ruido, "model"].tolist()

# Guardar listas
pd.Series(model_ruido).to_csv(os.path.join(RESULTS_DIR, "lista_ruido.txt"),
                              index=False, header=False)
pd.Series(model_no_ruido).to_csv(os.path.join(RESULTS_DIR, "lista_no_ruido.txt"),
                                 index=False, header=False)
pd.DataFrame({"ruido": model_ruido}).to_csv(os.path.join(RESULTS_DIR, "lista_ruido.csv"), index=False)
pd.DataFrame({"no_ruido": model_no_ruido}).to_csv(os.path.join(RESULTS_DIR, "lista_no_ruido.csv"), index=False)

print(f"Total estructuras: {len(df_proj)}")
print(f"Ruido: {len(model_ruido)}  |  No ruido: {len(model_no_ruido)}")
print("Ejemplos ruido:", model_ruido[:8])
print("Ejemplos no ruido:", model_no_ruido[:8])

# --- 2) Rehacer PCA SOLO con las no-ruido ---
# Importante: df_proj está ordenado alfabéticamente, pero X/labels están en el orden original.
# Construimos un reordenamiento para alinear X con df_proj:
model_to_idx = {m: i for i, m in enumerate(labels)}    # labels viene de arriba
order_idx = [model_to_idx[m] for m in df_proj["model"]] # X en el orden de df_proj
X_sorted = X[order_idx, :]

# Seleccionamos solo no-ruido
X_clean = X_sorted[~mask_ruido.values, :]
models_clean = df_proj.loc[~mask_ruido, "model"].tolist()

if X_clean.shape[0] < 2:
    raise RuntimeError("Quedaron menos de 2 estructuras no-ruido: no se puede rehacer el PCA.")

# PCA 2D (sklearn si está; si no, SVD)
try:
    from sklearn.decomposition import PCA
    pca_clean = PCA(n_components=2)
    proj_clean = pca_clean.fit_transform(X_clean)
    var1_c, var2_c = pca_clean.explained_variance_ratio_[0], pca_clean.explained_variance_ratio_[1]
except Exception:
    Xc2 = X_clean - X_clean.mean(axis=0, keepdims=True)
    U2, S2, Vt2 = np.linalg.svd(Xc2, full_matrices=False)
    proj_clean = U2[:, :2] * S2[:2]
    var = (S2**2) / (X_clean.shape[0]-1)
    var_ratio2 = var / var.sum()
    var1_c, var2_c = float(var_ratio2[0]), float(var_ratio2[1])

df_clean = pd.DataFrame({"model": models_clean, "PC1": proj_clean[:,0], "PC2": proj_clean[:,1]})
df_clean.to_csv(os.path.join(RESULTS_DIR, "pca_projection_kabsch_NO_RUIDO.csv"), index=False)

# --- 3) Etiquetas (1AKE_A y 4AKE_A con estrella; el resto numeradas) ---
SPECIAL = {"1AKE_A", "4AKE_A"}
display_clean = []
num_map_clean = []
counter = 1
for name in df_clean["model"]:
    if name in SPECIAL:
        display_clean.append(name)
    else:
        display_clean.append(str(counter))
        num_map_clean.append((counter, name))
        counter += 1
df_clean["display"] = display_clean

if num_map_clean:
    pd.DataFrame(num_map_clean, columns=["label_number", "model"]).to_csv(
        os.path.join(RESULTS_DIR, "pca_label_map_NO_RUIDO.csv"), index=False
    )

# --- 4) Gráfico solo PC1 vs PC2 (no-ruido) ---
plt.figure(figsize=(7, 7))

is_special_clean = df_clean["model"].isin(SPECIAL)
plt.scatter(df_clean.loc[~is_special_clean, "PC1"], df_clean.loc[~is_special_clean, "PC2"], alpha=0.85)
plt.scatter(df_clean.loc[is_special_clean, "PC1"],  df_clean.loc[is_special_clean,  "PC2"],
            s=80, marker="*", zorder=3)

for _, r in df_clean.iterrows():
    txt = r["display"]
    kws = dict(xytext=(3, 3), textcoords="offset points", fontsize=8)
    if r["model"] in SPECIAL:
        kws.update(fontweight="bold")
    plt.annotate(txt, (r["PC1"], r["PC2"]), **kws)

plt.xlabel(f"PC1 (noise removed; {var1_c*100:.1f}% var)")
plt.ylabel(f"PC2 (noise removed; {var2_c*100:.1f}% var)")
plt.title("Cα PCA (non-noise structures only)")
plt.grid(True, alpha=0.3)
plt.tight_layout()
out_png_clean = os.path.join(RESULTS_DIR, "pca_projection_kabsch_NO_RUIDO_PC12.png")
plt.savefig(out_png_clean, dpi=150)
plt.show()

print("\nListas guardadas en:")
print(" - lista_ruido.txt / .csv")
print(" - lista_no_ruido.txt / .csv")
print("PCA no-ruido:")
print(f" - pca_projection_kabsch_NO_RUIDO.csv")
print(f" - pca_label_map_NO_RUIDO.csv (si hay numerados)")
print(f" - {out_png_clean}")
