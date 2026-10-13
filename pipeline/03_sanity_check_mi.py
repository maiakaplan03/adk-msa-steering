#!/usr/bin/env python3
"""
03_sanity_check_mi.py — "3.3 métrica de coevolución: MI" and "3.4 control con query y MSA sintéticos".

Builds a synthetic ("toy") alignment of 120 sequences x 40 columns in which columns 5 and
30 are coupled (70% of rows follow a fixed residue mapping), plus a small contact list with
a hub on residue 5. Then runs every action x mode over 40 seeds and checks that the edit
moves the mutual information of the edited pairs in the expected direction:
promote / promote_joint should raise MI, break should lower it.

Writes results/sanity_check/mi_delta_mean_sem_toy.png and mi_direction_by_action_mode.png.
Usage: python pipeline/03_sanity_check_mi.py   (runs in a scratch folder ./demo_B_multiseed)
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from msa_editor import write_fasta, read_fasta_or_a3m, make_progressive_msas_weighted  # noqa: E402

OUT_FIG_DIR = os.environ.get("MI_SANITY_OUT", os.path.join(os.path.dirname(HERE), "results", "sanity_check"))
os.makedirs(OUT_FIG_DIR, exist_ok=True)


def plot_sanity_direction_consistency(df, use_metric="delta_mi_mean_all", savepath=None):
    actions = ["break", "promote", "promote_joint"]
    modes   = ["uniform", "degree", "hybrid"]

    tmp = df.copy()

    def is_ok(row):
        x = row[use_metric]
        if pd.isna(x):
            return np.nan
        if row["action"] == "break":
            return x < 0
        else:
            return x > 0

    tmp["ok"] = tmp.apply(is_ok, axis=1)

    frac_ok = (
        tmp.pivot_table(index="action", columns="mode", values="ok", aggfunc="mean")
           .reindex(actions)[modes]
    )

    plt.figure(figsize=(7.2, 3.6))
    ax = plt.gca()
    im = ax.imshow(frac_ok.values, aspect="auto", vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(modes))); ax.set_xticklabels(modes)
    ax.set_yticks(np.arange(len(actions))); ax.set_yticklabels(actions)
    ax.set_xlabel("MODE"); ax.set_ylabel("ACT")
    ax.set_title("Sanity check — fraction of seeds with expected ΔMI direction")

    # SOLO flechas
    for i, act in enumerate(actions):
        for j, mode in enumerate(modes):
            v = frac_ok.values[i, j]
            if np.isnan(v):
                txt = "–"
            else:
                txt = "↓" if act == "break" else "↑"
            ax.text(j, i, txt, ha="center", va="center", fontsize=18)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("fraction OK (0–1)")
    plt.tight_layout()
    if savepath:
        plt.savefig(savepath, dpi=300, bbox_inches="tight")
    plt.show()



# ---------- 1) Toy moderado (promote puede mejorar) ----------
def make_toy_msa_and_contacts_moderate(out_dir="demo_toyB", nseq=120, L=40, seed=0, p_corr=0.70):
    rng = np.random.default_rng(seed)
    os.makedirs(out_dir, exist_ok=True)

    aa = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    target_pair = (5, 30)  # 1-based

    # Mapeo “verdadero”
    map_AB = {"K": "E", "R": "D", "A": "A", "V": "V", "G": "G"}
    keys = np.array(list(map_AB.keys()))
    vals = np.array(list(map_AB.values()))

    headers, seqs = [], []
    for s in range(nseq):
        x = rng.choice(aa, size=L)

        a = rng.choice(keys)
        x[target_pair[0]-1] = a

        # 70% correlacionado, 30% “ruido” (random)
        if rng.random() < p_corr:
            x[target_pair[1]-1] = map_AB[str(a)]
        else:
            x[target_pair[1]-1] = rng.choice(aa)

        headers.append(f"toy_{s}")
        seqs.append("".join(x.tolist()))

    msa_path = os.path.join(out_dir, "toy.a3m")
    write_fasta(msa_path, headers, seqs)

    # CSV de contactos: hacemos hub en 5 + extras
    hub_i = target_pair[0]
    hub_partners = [12, 18, 22, 26, 30, 33, 35]
    extra_pairs = [(8, 20), (10, 27), (14, 37), (16, 25), (19, 36), (2, 39)]
    pairs = [(hub_i, j) for j in hub_partners] + extra_pairs

    # distancias para que "hybrid" exista
    dist = []
    for (i, j) in pairs:
        if (i, j) == target_pair or (j, i) == target_pair:
            dist.append(6.0)
        elif i == hub_i:
            dist.append(float(rng.uniform(7.0, 12.0)))
        else:
            dist.append(float(rng.uniform(9.0, 15.0)))

    contacts_path = os.path.join(out_dir, "toy_contacts.csv")
    pd.DataFrame({"resi1":[p[0] for p in pairs], "resi2":[p[1] for p in pairs], "dist_B":dist}).to_csv(contacts_path, index=False)

    ref_index = 0
    ref_seq = seqs[ref_index]
    return msa_path, contacts_path, ref_index, target_pair, ref_seq


# ---------- helpers lectura de resultados ----------
def _read_one_run(out_dir):
    summ = pd.read_csv(os.path.join(out_dir, "summary.csv"))
    tag = summ.loc[0, "tag"]
    rep = pd.read_csv(os.path.join(out_dir, f"report_{tag}.csv"))
    return tag, summ, rep

def _extract_pair_delta(rep_df, target_pair):
    ti, tj = target_pair
    tgt = (min(ti, tj), max(ti, tj))
    for _, row in rep_df.iterrows():
        s = str(row["pair"]).strip().strip("()")
        parts = [x.strip() for x in s.split(",")]
        if len(parts) != 2:
            continue
        try:
            a = int(parts[0]); b = int(parts[1])
        except:
            continue
        cur = (min(a, b), max(a, b))
        if cur == tgt:
            return float(row["delta_mi"]), float(row["mi_before"]), float(row["mi_after"])
    return None


# ---------- 2) Multi-seed experiment ----------
def run_multi_seed_demo_B(
    nseq=120, L=40, p_corr=0.70,
    pair_frac=0.20, strength=0.60,
    alpha=0.05, gamma=4.0,
    seeds=range(1, 41),          # 40 replicados (ajustá a gusto)
):
    demo_root = "demo_B_multiseed"
    toy_dir = os.path.join(demo_root, "toy_data")
    out_root = os.path.join(demo_root, "runs")
    os.makedirs(out_root, exist_ok=True)

    msa_path, contacts_csv, ref_index, target_pair, ref_before = make_toy_msa_and_contacts_moderate(
        out_dir=toy_dir, nseq=nseq, L=L, seed=0, p_corr=p_corr
    )

    actions = ["break", "promote", "promote_joint"]
    modes   = ["uniform", "degree", "hybrid"]

    records = []
    for action in actions:
        for mode in modes:
            for sd in seeds:
                run_out_dir = os.path.join(out_root, f"act-{action}__mode-{mode}__seed-{sd}")

                out_dir = make_progressive_msas_weighted(
                    msa_in_path=msa_path,
                    contacts_csv=contacts_csv,
                    base_out_dir=run_out_dir,
                    mode=mode,
                    pair_fracs=[pair_frac],
                    strength_grid=[strength],
                    seeds=[sd],
                    action=action,
                    ref_index=ref_index,
                    one_based=True,
                    alpha=alpha,
                    power_gamma=gamma,
                    min_seq_sep=0,
                    use_dist_filter=False,
                    max_dist_for_gain=None
                )

                tag, summ, rep = _read_one_run(out_dir)
                # % de pares editados con ΔMI > 0 en esta corrida
                deltas = rep["delta_mi"].astype(float).values
                pct_pos = float(np.mean(deltas > 0)) if len(deltas) else np.nan


                # mean ΔMI sobre pares elegidos
                dmi_mean = float(summ.loc[0, "delta_mi_mean"]) if not pd.isna(summ.loc[0, "delta_mi_mean"]) else np.nan

                # ref intacta (sanity)
                msa_out_path = os.path.join(out_dir, f"msa_{tag}.a3m")
                _, seqs_out = read_fasta_or_a3m(msa_out_path)
                ref_ok = (seqs_out[ref_index] == ref_before)

                records.append({
                    "action": action, "mode": mode, "seed": sd,
                    "delta_mi_mean_all": dmi_mean,
                    "pct_pos": pct_pos,
                    "ref_intact": ref_ok
                })

    df = pd.DataFrame(records)
    df.to_csv(os.path.join(demo_root, "all_results.csv"), index=False)
    return df, demo_root


# ---------- 3) Figuras 1 y 3 con media ± SEM ----------
def sem(x):
    x = np.asarray(x, dtype=float)
    x = x[~np.isnan(x)]
    if len(x) <= 1:
        return np.nan
    return x.std(ddof=1) / np.sqrt(len(x))

def plot_fig1_heatmap(df):
    actions = ["break", "promote", "promote_joint"]
    modes   = ["uniform", "degree", "hybrid"]

    mean_piv = df.pivot_table(index="action", columns="mode", values="delta_mi_mean_all", aggfunc="mean").reindex(actions)[modes]
    sem_piv  = df.pivot_table(index="action", columns="mode", values="delta_mi_mean_all", aggfunc=sem).reindex(actions)[modes]

    fig = plt.figure(figsize=(7.8, 4.0))
    ax = plt.gca()
    im = ax.imshow(mean_piv.values, aspect="auto")

    ax.set_xticks(np.arange(len(modes))); ax.set_xticklabels(modes)
    ax.set_yticks(np.arange(len(actions))); ax.set_yticklabels(actions)
    ax.set_title("Figure 1 — mean(ΔMI over edited pairs) ± SEM (multi-seed)")
    ax.set_xlabel("MODE"); ax.set_ylabel("ACT")

    for i in range(len(actions)):
        for j in range(len(modes)):
            m = mean_piv.values[i, j]
            s = sem_piv.values[i, j]
            ax.text(j, i, f"{m:.3f}±{s:.3f}", ha="center", va="center")

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("mean ΔMI")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_FIG_DIR, "mi_delta_mean_sem_toy.png"), dpi=200, bbox_inches="tight")
    plt.show()


# ---------- RUN ----------
dfB, root = run_multi_seed_demo_B(
    nseq=120, L=40, p_corr=0.70,
    pair_frac=0.20, strength=0.60,
    alpha=0.05, gamma=4.0,
    seeds=range(1, 41)   # 40 seeds
)

print("Sanity check: ref_intact all True?", dfB["ref_intact"].all())
plot_fig1_heatmap(dfB)
plot_sanity_direction_consistency(dfB, savepath=os.path.join(OUT_FIG_DIR, "mi_direction_by_action_mode.png"))
print("Resultados guardados en:", os.path.abspath(root))
