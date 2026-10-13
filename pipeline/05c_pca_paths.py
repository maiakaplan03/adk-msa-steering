#!/usr/bin/env python3
"""
05c_pca_paths.py — paths through the PCA plane (MST over a k-NN graph) and their
simplification to ~15 representative structures per path.

Usage:
  python pipeline/05c_pca_paths.py --pca-csv results_pca/pca_projection_kabsch_NO_RUIDO.csv \
         --map-csv results_pca/pca_label_map_NO_RUIDO.csv --out-dir pca_paths_results \
         [--end2 referencia]
"""
import argparse
_ap = argparse.ArgumentParser()
_ap.add_argument("--pca-csv", default="pca_projection_kabsch_NO_RUIDO.csv")
_ap.add_argument("--map-csv", default="pca_label_map_NO_RUIDO.csv")
_ap.add_argument("--out-dir", default="pca_paths_results")
_ap.add_argument("--start", default="1AKE_A", help="model name of the start anchor")
_ap.add_argument("--end1", default="4AKE_A", help="model name of the first end anchor")
_ap.add_argument("--end2", default="referencia", help="model name of the second end point")
_ap.add_argument("--k", type=int, default=18)
_ap.add_argument("--max-points", type=int, default=15)
_args = _ap.parse_args()
import matplotlib
matplotlib.use("Agg")

# ============================================================
# 🧭 Caminos en PCA - SOLO MST(kNN)
# - Lee:
#   /content/pca_projection_kabsch_NO_RUIDO.csv  -> model, PC1, PC2
#   /content/pca_label_map_NO_RUIDO.csv          -> label_number, model
#
# Genera 2 caminos:
#   1) del punto más cercano a 1AKE_A al punto más cercano a 4AKE_A
#   2) del punto más cercano a 1AKE_A al punto más cercano a 'referencia'
# ============================================================
# Requisitos (una vez):
# !pip -q install networkx

import os, math, numpy as np, pandas as pd, matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import networkx as nx

# ---------- RUTAS (NO Drive) ----------
PC_CSV  = _args.pca_csv
MAP_CSV = _args.map_csv
OUT_DIR = _args.out_dir
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- CONFIG: puntos de referencia ----------
# Estos son los "anclajes" (cristalinos / puntos especiales)
ANCHOR_START_NAME = _args.start
ANCHOR_END1_NAME  = _args.end1
ANCHOR_END2_NAME  = _args.end2  # debe existir en la columna 'model'

# Si quisieras usarlos por coordenadas en lugar de nombres, poné None arriba
# y definí aquí las coordenadas PCA (PC1, PC2):
ANCHOR_START_COORD = None  # ej. (62.8, -31.7)
ANCHOR_END1_COORD  = None
ANCHOR_END2_COORD  = None

# ---------- Parámetros de grafo ----------
K_NEIGHBORS = _args.k     # 15–25 suele ir bien
USE_MUTUAL  = False  # True = kNN mutuo (más robusto, menos aristas)

# ---------- CARGA ----------
df = pd.read_csv(PC_CSV)          # model, PC1, PC2
map_df = pd.read_csv(MAP_CSV)     # label_number, model

# unir para tener label_number junto con PC1/PC2
df = df.merge(map_df[['label_number','model']], on='model', how='left')

# matriz de puntos y nombres
X = df[['PC1','PC2']].to_numpy(float)
names = df['model'].astype(str).to_numpy()

# ---------- Helpers para índices/coordenadas ----------
def index_by_name(name):
    hits = np.where(names == str(name))[0]
    return int(hits[0]) if len(hits) else None

def coord_from_name(name):
    idx = index_by_name(name)
    if idx is None:
        raise ValueError(f"No se encontró el nombre '{name}' en la columna 'model'.")
    return X[idx]

def nearest_idx_to_anchor(anchor_name=None, anchor_coord=None, exclude_name=None):
    """
    Devuelve el índice del punto más cercano al ancla:
    - Si anchor_name está definido: toma su coordenada en PCA.
    - Si no, usa anchor_coord (tupla PC1,PC2).
    - Si exclude_name no es None, excluye todos los puntos con ese nombre
      (útil para que el "punto más cercano a 1AKE" NO sea 1AKE en sí).
    """
    if anchor_name is not None:
        base_xy = coord_from_name(anchor_name)
    elif anchor_coord is not None:
        base_xy = np.asarray(anchor_coord, float)
    else:
        raise ValueError("Definí anchor_name o anchor_coord.")

    base_xy = np.asarray(base_xy, float)
    mask = np.ones(len(X), dtype=bool)
    if exclude_name is not None:
        mask &= (names != str(exclude_name))

    if not np.any(mask):
        raise RuntimeError("La máscara de exclusión dejó 0 puntos.")

    X_masked = X[mask]
    d2 = np.sum((X_masked - base_xy)**2, axis=1)
    rel_idx = int(np.argmin(d2))
    global_indices = np.arange(len(X))[mask]
    return int(global_indices[rel_idx])

# ---------- Grafo k-NN ----------
def build_knn_edges(X, k=15, mutual=False):
    n = len(X)
    k = min(k+1, n)  # incluye self en la col 0
    nn = NearestNeighbors(n_neighbors=k).fit(X)
    dists, idxs = nn.kneighbors(X)
    edges = []
    for i in range(n):
        for d, j in zip(dists[i,1:], idxs[i,1:]):  # omitir self
            if mutual and i not in idxs[j,1:]:
                continue
            edges.append((int(i), int(j), float(d)))
    return np.array(edges, dtype=float), idxs, dists

edges, idxs, dists = build_knn_edges(X, K_NEIGHBORS, USE_MUTUAL)

# ---------- Método MST(kNN) ----------
def path_mst_knn(X, edges, s, t):
    G = nx.Graph()
    for i, j, d in edges:
        G.add_edge(int(i), int(j), weight=float(d))
    T = nx.minimum_spanning_tree(G, weight="weight")
    if not nx.has_path(T, s, t):
        raise RuntimeError("Sin camino en el MST (subí K_NEIGHBORS).")
    return nx.shortest_path(T, s, t, weight="weight")

# ---------- Resolver índices de los puntos "más cercanos" ----------
# Punto más cercano a 1AKE_A (excluyendo el propio 1AKE_A)
start_idx = nearest_idx_to_anchor(
    anchor_name=ANCHOR_START_NAME,
    anchor_coord=ANCHOR_START_COORD,
    exclude_name=ANCHOR_START_NAME
)

# Punto más cercano a 4AKE_A (excluyendo el propio 4AKE_A)
end1_idx = nearest_idx_to_anchor(
    anchor_name=ANCHOR_END1_NAME,
    anchor_coord=ANCHOR_END1_COORD,
    exclude_name=ANCHOR_END1_NAME
)

# Punto más cercano a 'referencia' (excluyendo el propio 'referencia')
end2_idx = nearest_idx_to_anchor(
    anchor_name=ANCHOR_END2_NAME,
    anchor_coord=ANCHOR_END2_COORD,
    exclude_name=ANCHOR_END2_NAME
)

print("Puntos seleccionados para los caminos:")
print(f" - cercano a {ANCHOR_START_NAME}: idx={start_idx}, name={names[start_idx]}")
print(f" - cercano a {ANCHOR_END1_NAME}:  idx={end1_idx},  name={names[end1_idx]}")
print(f" - cercano a {ANCHOR_END2_NAME}:  idx={end2_idx},  name={names[end2_idx]}")

# ---------- Calcular SOLO 2 caminos (MST) ----------
path_1AKE_to_4AKE = path_mst_knn(X, edges, start_idx, end1_idx)
path_1AKE_to_ref  = path_mst_knn(X, edges, start_idx, end2_idx)

def names_from_path(path_idx):
    return [names[i] for i in path_idx]

paths = {
    "path_1AKE_to_4AKE": path_1AKE_to_4AKE,
    "path_1AKE_to_ref":  path_1AKE_to_ref,
}
path_names = {k: names_from_path(v) for k, v in paths.items()}

# ---------- Guardar listas ----------
for tag, lst in path_names.items():
    pd.DataFrame({"order": range(len(lst)), "name": lst}).to_csv(
        f"{OUT_DIR}/{tag}_path_names.csv", index=False
    )

# ---------- Plot ----------
plt.figure(figsize=(8,7))
plt.scatter(X[:,0], X[:,1], s=12, alpha=0.5, label="predictions")

def draw_path(path, style, label):
    XY = X[np.array(path)]
    plt.plot(XY[:,0], XY[:,1], style, linewidth=2.0, label=label)
    plt.scatter(XY[:,0], XY[:,1], s=22)

draw_path(path_1AKE_to_4AKE, '-',  "1AKE→4AKE (MST kNN)")
draw_path(path_1AKE_to_ref,  '--', "1AKE→reference (MST kNN)")

# marcamos los tres anclajes originales (cristales / puntos especiales)
for name, marker, color in [
    (ANCHOR_START_NAME, 'o', 'black'),
    (ANCHOR_END1_NAME,  's', 'black'),
    (ANCHOR_END2_NAME,  'D', 'black')
]:
    idx = index_by_name(name)
    if idx is not None:
        plt.scatter([X[idx,0]], [X[idx,1]], s=120, marker=marker, label=f"anchor: {name}")

plt.title("Paths in the PCA plane (kNN MST only)")
plt.xlabel("PC1"); plt.ylabel("PC2")
plt.legend(loc="best", fontsize=9)
plt.grid(True, alpha=0.3)
fig_path = f"{OUT_DIR}/pca_paths_overlay.png"
plt.savefig(fig_path, dpi=200, bbox_inches="tight")
plt.show()

print("✅ Resultados en:", OUT_DIR)
for k, v in paths.items():
    print(f" - {k}: {len(v)} nodos -> {OUT_DIR}/{k}_path_names.csv")
print("✅ Figura:", fig_path)


# ============================================================
# 🧭 Post-proceso de caminos MST: 15 estructuras sin cristales
# Usa como INPUT los outputs del bloque anterior:
#   /content/pca_paths_2_results/path_1AKE_to_4AKE_path_names.csv
#   /content/pca_paths_2_results/path_1AKE_to_ref_path_names.csv
#
# Además necesita la proyección PCA para obtener coordenadas:
#   /content/pca_projection_kabsch_NO_RUIDO(2).csv  (model, PC1, PC2)
#
# Condiciones:
#   - 15 estructuras representativas por camino (≈ MAX_POINTS)
#   - Excluir 1AKE_A y 4AKE_A de los grupos finales
#   - Conservar las dos puntas del camino (inicio y fin, ya sin cristales)
#   - Mostrar gráficamente qué estructuras quedaron en cada grupo
# ============================================================

import os, numpy as np, pandas as pd, matplotlib.pyplot as plt

# ---------- Rutas (deben coincidir con el bloque anterior) ----------
PC_CSV   = _args.pca_csv
OUT_DIR  = _args.out_dir
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- Parámetros ----------
MAX_POINTS     = _args.max_points
BANNED_NAMES   = {"1AKE_A", "4AKE_A"}   # cristales a excluir de los 15 finales

# ---------- Carga de datos ----------
df = pd.read_csv(PC_CSV)  # columnas: model, PC1, PC2
# Mapa nombre->(PC1,PC2)
coord = df.set_index("model")[["PC1","PC2"]].to_dict("index")

def load_path_names(tag):
    """Lee un camino crudo generado antes, devuelve lista de nombres."""
    path_csv = f"{OUT_DIR}/{tag}_path_names.csv"
    if not os.path.exists(path_csv):
        raise FileNotFoundError(f"No encuentro {path_csv}. ¿Corriste el bloque anterior?")
    return pd.read_csv(path_csv)["name"].astype(str).tolist()

def names_to_xy(names):
    """Convierte lista de modelos a array Nx2 de coordenadas (PC1,PC2)."""
    pts = []
    missing = []
    for n in names:
        if n in coord:
            pts.append([coord[n]["PC1"], coord[n]["PC2"]])
        else:
            missing.append(n)
    if missing:
        print(f"⚠️ Aviso: {len(missing)} nombres no tienen coordenadas y se omiten (ej: {missing[:3]} ...)")
    return np.array(pts, float)

def simplify_uniform_keep_endpoints(names, max_points):
    """
    Selección aproximadamente uniforme a lo largo del camino,
    CONSERVANDO SIEMPRE:
      - el primer elemento (punto inicial)
      - el último elemento (punto final)
    """
    if (max_points is None) or (len(names) <= max_points):
        return names

    XY = names_to_xy(names)
    if len(XY) <= 2:
        return names

    # longitud acumulada a lo largo del camino en el espacio PC1/PC2
    seg = np.sqrt(np.sum(np.diff(XY, axis=0)**2, axis=1))
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    total = cum[-1]
    targets = np.linspace(0, total, max_points)

    out = []
    for t in targets:
        i = np.searchsorted(cum, t, side="left")
        i = min(i, len(names)-1)
        out.append(names[i])

    # asegurar extremos y quitar consecutivos duplicados
    out[0]  = names[0]
    out[-1] = names[-1]
    uniq = [out[0]]
    for v in out[1:]:
        if v != uniq[-1]:
            uniq.append(v)
    return uniq

def path_length(names):
    XY = names_to_xy(names)
    if len(XY) < 2:
        return 0.0
    return float(np.sum(np.sqrt(np.sum(np.diff(XY,axis=0)**2,axis=1))))

# ---------- Procesa SOLO LOS 2 CAMINOS ----------
tags = ["path_1AKE_to_4AKE", "path_1AKE_to_ref"]

paths_raw = {}
paths_no_crystals = {}
paths_simpl = {}

for t in tags:
    raw = load_path_names(t)
    paths_raw[t] = raw

    # 1) quitamos cristales (pero mantenemos orden)
    no_crystals = [n for n in raw if n not in BANNED_NAMES]
    if len(no_crystals) < 2:
        raise RuntimeError(f"El camino {t} se quedó sin puntos suficientes tras quitar cristales.")
    paths_no_crystals[t] = no_crystals

    # 2) simplificación uniforme sobre el camino SIN cristales
    simpl = simplify_uniform_keep_endpoints(no_crystals, MAX_POINTS)
    paths_simpl[t] = simpl

    # chequeo: puntas del camino sin cristales deben estar en el grupo simplificado
    if no_crystals[0] not in simpl or no_crystals[-1] not in simpl:
        raise RuntimeError(f"La simplificación del camino {t} perdió una de las puntas (sin cristales).")

    # chequeo: que los cristales no estén en el grupo final
    bad = [n for n in simpl if n in BANNED_NAMES]
    if bad:
        raise RuntimeError(f"El camino simplificado {t} aún contiene cristales: {bad}")

# ---------- Guarda CSV de caminos simplificados ----------
for t in tags:
    out_csv = f"{OUT_DIR}/{t}_path_names_simplified.csv"
    pd.DataFrame({
        "order": range(len(paths_simpl[t])),
        "name":  paths_simpl[t]
    }).to_csv(out_csv, index=False)

# ---------- Plot: muestra qué estructuras quedaron en cada grupo ----------
plt.figure(figsize=(8,7))
plt.scatter(df["PC1"], df["PC2"], s=10, alpha=0.25, label="predictions")

def draw_by_names(names, style, label, lw=2.0, scatter_size=25):
    XY = names_to_xy(names)
    if len(XY) >= 2:
        plt.plot(XY[:,0], XY[:,1], style, linewidth=lw, label=label)
    plt.scatter(XY[:,0], XY[:,1], s=scatter_size)

# colores/estilos para distinguir caminos y simplificados
# (no importa el color concreto, solo que sean distintos)
styles_raw = {
    "path_1AKE_to_4AKE": ":",
    "path_1AKE_to_ref":  "--",
}
styles_simpl = {
    "path_1AKE_to_4AKE": "-",
    "path_1AKE_to_ref":  "-",
}

# Caminos crudos (sin quitar cristales) – líneas finas
for t in tags:
    draw_by_names(paths_raw[t],
                  styles_raw[t],
                  f"{t} raw (n={len(paths_raw[t])})",
                  lw=1.2, scatter_size=0)  # no marcamos puntos crudos

# Caminos simplificados (ya sin cristales) – líneas gruesas + puntos grandes
for t in tags:
    draw_by_names(paths_simpl[t],
                  styles_simpl[t],
                  f"{t} simpl (n={len(paths_simpl[t])})",
                  lw=2.8, scatter_size=40)

plt.title("MST paths in PCA: groups of 15 structures (without 1AKE_A / 4AKE_A)")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.legend(loc="best", fontsize=9)
plt.grid(True, alpha=0.3)

fig_path = f"{OUT_DIR}/pca_paths_simplified_groups.png"
plt.savefig(fig_path, dpi=220, bbox_inches="tight")
plt.show()

# ---------- Resumen ----------
print("✅ Resultados en:", OUT_DIR)
for t in tags:
    Lr = path_length(paths_no_crystals[t])
    Ls = path_length(paths_simpl[t])
    print(f" - {t:16s}: sin cristales n={len(paths_no_crystals[t]):4d}, len={Lr:7.2f} | "
          f"simpl n={len(paths_simpl[t]):4d}, len={Ls:7.2f} -> "
          f"{OUT_DIR}/{t}_path_names_simplified.csv")

print("✅ Figura:", fig_path)
print(f"(Simplificación) MAX_POINTS={MAX_POINTS}, excluyendo {BANNED_NAMES}")
print("\nEjemplo nombres seleccionados (primer camino):")
print("path_1AKE_to_4AKE:", paths_simpl["path_1AKE_to_4AKE"])
