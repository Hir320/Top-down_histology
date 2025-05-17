# =============================================================
# ACC–V1 Topography & Bias Analysis Pipeline
# -------------------------------------------------------------
# Author : ChatGPT (refactor for Oishi lab)
# Date   : 2025‑05‑16
# Purpose: Re‑implement Morimoto et al. (2021) style topography/
#          bias analysis for tangential 10× volumes without CCF.
# =============================================================

"""
Usage (Jupyter Lab):

```python
%run acc_v1_topography_analysis.py  # <- define utils

# --- parameters --------------------------------------------------
animals  = [1]                 # CTB IDs as integers
regions  = ["ACC"]            # add e.g. "RSC" if available
root4x   = Path("./4X")
root10x  = Path("./Imaris")

# --- analysis ----------------------------------------------------
all_stats = []
for aid in animals:
    for reg in regions:
        stats, figs = analyze_one(animal_id=aid,
                                  region=reg,
                                  root4x=root4x,
                                  root10x=root10x,
                                  plot=True)
        all_stats.append(stats)
```

Functions will emit Plotly figures inline and return a `stats` dict
containing Procrustes 残差, 1‑D ρ, bias χ², permutation p values, etc.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict

from sklearn.decomposition import PCA
from sklearn.utils import shuffle
from scipy.spatial import procrustes
from scipy.stats import pearsonr, chi2

import plotly.graph_objects as go
import plotly.express as px

# -----------------------------------------------------------------
# 0.  Inject‑site → 共通座標系変換 (already refined)
# -----------------------------------------------------------------

EXPECT_RIGHTWARD = {"CTB1": False, "CTB2": True, "CTB3": True, "CTB4": False}
ML_SIGN          = {"CTB1": -1}


def get_transform(animal: str, root4x: Path = Path("./4X")) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (v_AP[2], v_ML[2], origin[3], R[3×3]). Units = µm."""
    labels = (pd.read_csv(root4x / "inj_label.csv")
                .set_index(["animal", "channel"]))
    inj_df = (pd.read_csv(root4x / f"{animal}_4X_1.csv")
                .assign(channel=lambda d: d["Label"].str.extract(r"C(\d+)").astype(int))
                .merge(labels.loc[animal].reset_index(), on="channel"))

    order_ap = {"ant": 0, "mid": 1, "pos": 2, "post": 2}
    inj_df   = inj_df.assign(_ap=inj_df.AP.map(order_ap)).sort_values("_ap")

    xy       = inj_df[["X", "Y"]].to_numpy(float)
    ant_xy   = inj_df.loc[inj_df.AP == "ant",  ["X", "Y"]].to_numpy(float)
    pos_xy   = inj_df.loc[inj_df.AP.isin(["pos", "post"]), ["X", "Y"]].to_numpy(float)

    v_AP = (pos_xy - ant_xy).ravel();  v_AP /= np.linalg.norm(v_AP)
    v_ML = np.array([-v_AP[1], v_AP[0]])

    if ((inj_df.loc[inj_df.ML == "lat", ["X", "Y"]].to_numpy(float) - xy.mean(0)) @ v_ML).mean() < 0:
        v_ML *= -1

    if EXPECT_RIGHTWARD.get(animal, True) != (v_AP[0] > 0):
        v_AP *= -1  # v_ML stay

    v_ML *= ML_SIGN.get(animal, 1)

    origin = np.append(xy.mean(0), 0.0)

    R = np.eye(3, dtype=float)
    R[:2, 0] = v_AP
    R[:2, 1] = v_ML
    return v_AP, v_ML, origin, R

# -----------------------------------------------------------------
# 1.  I/O helpers
# -----------------------------------------------------------------

def load_injection(animal_id: int, root4x: Path) -> np.ndarray:
    """Return 3×2 array in µm: ant, mid, post injection centroids."""
    csv = root4x / f"CTB{animal_id}_4X_1.csv"
    df  = pd.read_csv(csv)
    for xc, yc in (("XM", "YM"), ("X", "Y")):
        if xc in df.columns and yc in df.columns:
            arr = df[[xc, yc]].astype(float).values
            break
    else:
        raise KeyError("Injection CSV must have X/Y or XM/YM")
    return arr


# def load_region_cells(animal_id: int, region: str, root10x: Path) -> np.ndarray:
#     """Concatenate all cells (all three channels) in given region. Return N×3."""
#     pts = []
#     for ch in (1, 2, 3):
#         f = root10x / f"CTB{animal_id}_{region}_10X_Ch{ch}_Position.csv"
#         if not f.exists():
#             continue
#         df = pd.read_csv(f, skiprows=3)
#         pts.append(df[["Position X", "Position Y", "Position Z"]].to_numpy(float))
#     if not pts:
#         raise FileNotFoundError("No region files found")
#     return np.vstack(pts)

def load_region_cells(animal_id: int, region: str, root10x: Path):
    """各チャネルを区別して読み込み → (N,3) xyz と (N,) labels"""
    xyz_all, lbl_all = [], []
    for ch, lbl in zip((1, 2, 3), (0, 1, 2)):          # Ch1→0, Ch2→1, Ch3→2
        f = root10x / f"{animal_id}_{region}_10X_Ch{ch}_Position.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f, skiprows=3)
        xyz_all.append(df[["Position X", "Position Y", "Position Z"]].to_numpy(float))
        lbl_all.append(np.full(len(df), lbl))
    xyz = np.vstack(xyz_all)
    labels = np.concatenate(lbl_all).astype(int)
    return xyz, labels        # ← ラベル付きで返す


# -----------------------------------------------------------------
# 2.  Geometry helpers
# -----------------------------------------------------------------

def pca_flatten(xyz: np.ndarray) -> Tuple[np.ndarray, PCA]:
    """Return N×2 flattened coords and fitted PCA object."""
    pca = PCA(n_components=2).fit(xyz)
    uv  = pca.transform(xyz)
    return uv, pca


def assign_color_labels(acc_uv: np.ndarray, inj_uv: np.ndarray) -> np.ndarray:
    """Return label idx (0/1/2) for each ACC cell via nearest injection in uv‑space."""
    dists = ((acc_uv[:, None, :] - inj_uv[None, :, :]) ** 2).sum(-1)  # N×3
    return dists.argmin(1)

# -----------------------------------------------------------------
# 3.  Metrics
# -----------------------------------------------------------------

def procrustes_disp(inj: np.ndarray, proj: np.ndarray) -> float:
    _, _, disp = procrustes(inj, proj)
    return disp


def corr_1d(inj: np.ndarray, proj: np.ndarray) -> float:
    """Project both sets onto first SVD axis and compute Pearson r."""
    u = inj - inj.mean(0)
    v = proj - proj.mean(0)
    axis = np.linalg.svd(u)[2][0]
    inj_1d = u @ axis
    proj_1d = v @ axis
    r, _ = pearsonr(inj_1d, proj_1d)
    return r


def bias_chi2(counts: np.ndarray) -> Tuple[float, float]:
    """Chi‑square test for uniformity across 3 injection colors."""
    expected = counts.mean() * np.ones_like(counts)
    chi2_stat = ((counts - expected) ** 2 / expected).sum()
    p = 1 - chi2.cdf(chi2_stat, df=len(counts) - 1)
    return chi2_stat, p


# -----------------------------------------------------------------
# 4.  Permutation tests
# -----------------------------------------------------------------

def permutation_p(obs: float, acc_uv: np.ndarray, inj_uv: np.ndarray, n_iter: int = 2000) -> float:
    """Two‑sided permutation p for Procrustes disp or 1‑D corr."""
    labels = assign_color_labels(acc_uv, inj_uv)
    acc_uv_sh = acc_uv.copy()
    surrogates = []
    for _ in range(n_iter):
        np.random.shuffle(labels)  # in‑place shuffle
        centroids = np.vstack([acc_uv_sh[labels == k].mean(0) for k in range(3)])
        surrogates.append(procrustes_disp(inj_uv, centroids))
    surrogates = np.array(surrogates)
    return (np.sum(surrogates <= obs) + np.sum(surrogates >= obs)) / n_iter

# -----------------------------------------------------------------
# 5.  Visualization helpers (Plotly)
# -----------------------------------------------------------------

def plot_overlay(inj_uv: np.ndarray, acc_uv: np.ndarray, labels: np.ndarray, title: str) -> go.Figure:
    col = ["orange", "limegreen", "dodgerblue"]
    sym_inj = "triangle-up"; sym_proj = "x"
    fig = go.Figure()
    for k in range(3):
        fig.add_trace(go.Scatter(x=[inj_uv[k, 0]], y=[inj_uv[k, 1]],
                                 mode="markers+text", marker=dict(color=col[k], symbol=sym_inj, size=14),
                                 text=[f"inj-{k}"], textposition="top center",
                                 name=f"Injection {k}"))
        sel = acc_uv[labels == k]
        fig.add_trace(go.Scatter(x=sel[:, 0], y=sel[:, 1], mode="markers",
                                 marker=dict(color=col[k], symbol=sym_proj, size=6, opacity=0.4),
                                 name=f"Proj {k}"))
    fig.update_layout(title=title, xaxis=dict(scaleanchor="y"), yaxis=dict(), width=550, height=500)
    return fig

# -----------------------------------------------------------------
# 6.  High‑level wrapper
# -----------------------------------------------------------------

def analyze_one(animal_id: int,
                region: str,
                root4x: Path = Path("./4X"),
                root10x: Path = Path("./Imaris"),
                plot: bool = True) -> Tuple[Dict, List[go.Figure]]:
    """Run full pipeline for a single animal & region."""
    animal = f"CTB{animal_id}"
    # 6‑1 load + transform --------------------------------------------------
    v_AP, v_ML, origin, R = get_transform(animal, root4x)
    inj_xy = load_injection(animal_id, root4x)          # 3×2
    xyz_acc = load_region_cells(animal_id, region, root10x)

    acc_flat = (xyz_acc - origin) @ R[:, :2]            # N×2 µm

    # 6‑2 PCA --------------------------------------------------------------
    inj_uv, pca = pca_flatten(inj_xy)
    acc_uv      = pca.transform(acc_flat)

    # 6‑3 label assignment -----------------------------------------------
    labels = assign_color_labels(acc_uv, inj_uv)

    # 6‑4 metrics ----------------------------------------------------------
    proj_cent = np.vstack([acc_uv[labels == k].mean(0) for k in range(3)])
    disp = procrustes_disp(inj_uv, proj_cent)
    r1d = corr_1d(inj_uv, proj_cent)

    counts = np.bincount(labels, minlength=3)
    chi2_stat, p_chi = bias_chi2(counts)

    # permutation (optional; can be slow)
    # p_perm = permutation_p(disp, acc_uv, inj_uv, n_iter=1000)
    p_perm = np.nan

    # 6‑5 visualize --------------------------------------------------------
    figs = []
    if plot:
        figs.append(plot_overlay(inj_uv, acc_uv, labels,
                                 title=f"{animal} {region}<br>disp={disp:.3f}, r1d={r1d:.2f}"))
        # Bias bar
        fig_bar = px.bar(x=["ant", "mid", "post"], y=counts,
                         title="Projection counts per injection color",
                         labels=dict(x="Injection", y="Cell count"))
        figs.append(fig_bar)

    stats = dict(animal=animal, region=region, disp=disp, r1d=r1d,
                 chi2=chi2_stat, p_chi=p_chi, p_perm=p_perm,
                 n_cells=len(acc_uv))
    return stats, figs

# -----------------------------------------------------------------
if __name__ == "__main__":
    import argparse, json, webbrowser, tempfile
    parser = argparse.ArgumentParser(description="Run one‑off analysis on CLI")
    parser.add_argument("animal_id", type=int, help="e.g. 1 for CTB1")
    parser.add_argument("region", type=str, default="ACC")
    parser.add_argument("--root4x", default="./4X")
    parser.add_argument("--root10x", default="./Imaris")
    args = parser.parse_args()

    stats, figs = analyze_one(args.animal_id, args.region,
                              Path(args.root4x), Path(args.root10x), plot=True)
    print(json.dumps(stats, indent=2))

    # open plots in browser if run as script
    for fig in figs:
        html = tempfile.mktemp(suffix=".html")
        fig.write_html(html, auto_open=False)
        webbrowser.open(html)
