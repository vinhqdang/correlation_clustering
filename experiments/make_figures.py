"""The figures of paper_mpc/ that show the algorithms at work (support,
gap map, partition crossover).

    python experiments/make_figures.py            # all figures
    python experiments/make_figures.py crossover  # one of them

They are computed by running the solver code on a neighbourhood of ca-GrQc.
The result figures are written by experiments/mpc_tables.py."""
import glob
import json
import os
import sys
from collections import deque

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)
FIG = os.path.join(ROOT, "paper_mpc", "figures")
RES = os.path.join(ROOT, "results")

# categorical slots 1-3 of the reference palette (validated all-pairs), inks, grid
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
GRID, NEUTRAL = "#e5e4e0", "#c9c8c2"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.spines.top": False,
    "axes.spines.right": False, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.titlecolor": INK,
    "pdf.fonttype": 42, "font.family": "DejaVu Sans"})


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight")
    if os.environ.get("PREVIEW"):
        fig.savefig(os.path.join(os.environ["PREVIEW"], name + ".png"), dpi=130,
                    bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def grid(ax, axis="both"):
    ax.grid(True, axis=axis, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


# --------------------------------------------------------------------------
# the example instance: a BFS ball of 40 vertices in ca-GrQc
# --------------------------------------------------------------------------

def example(root=789, size=40):
    import datasets as D
    from ccbench.graph import induced_subgraph
    G = D.load("ca-GrQc")
    seen, q = {root}, deque([root])
    while q and len(seen) < size:
        u = q.popleft()
        for v in G.indices[G.indptr[u]:G.indptr[u + 1]]:
            if int(v) not in seen and len(seen) < size:
                seen.add(int(v))
                q.append(int(v))
    return induced_subgraph(G, np.array(sorted(seen)))


def layout(g):
    import networkx as nx
    H = nx.Graph()
    H.add_nodes_from(range(g.n))
    H.add_edges_from(map(tuple, g.edges()))
    pos = nx.kamada_kawai_layout(H)
    return np.array([pos[i] for i in range(g.n)])


def draw_graph(ax, g, P, lab=None, node_color=None, edge_color=None, size=34):
    """Edges inside a cluster dark, cut edges light; nodes coloured by node_color."""
    e = g.edges()
    for u, v in e:
        inside = lab is not None and lab[u] == lab[v]
        c = edge_color or (INK2 if inside else NEUTRAL)
        ax.plot(*P[[u, v]].T, color=c, lw=1.1 if (inside or edge_color) else 0.6,
                zorder=1, solid_capstyle="round", ls="-" if (inside or edge_color) else (0, (2, 1.5)))
    ax.scatter(P[:, 0], P[:, 1], s=size, c=node_color if node_color is not None else MUTED,
               edgecolors="white", linewidths=0.8, zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_aspect("equal")


def hulls(ax, P, lab, color="#dcdbd5", alpha=0.55):
    from scipy.spatial import ConvexHull
    for c in np.unique(lab):
        idx = np.flatnonzero(lab == c)
        if len(idx) < 2:
            continue
        pts = P[idx]
        if len(idx) >= 3:
            try:
                h = ConvexHull(pts)
                poly = pts[h.vertices]
            except Exception:
                poly = pts
        else:
            poly = pts
        cen = poly.mean(axis=0)
        d = poly - cen
        nrm = np.linalg.norm(d, axis=1, keepdims=True) + 1e-9
        poly = poly + 0.045 * d / nrm
        from matplotlib.patches import Polygon
        ax.add_patch(Polygon(poly, closed=True, fc=color, ec="none", alpha=alpha, zorder=0,
                             joinstyle="round"))


# --------------------------------------------------------------------------
# 1. distance-two support
# --------------------------------------------------------------------------

def fig_support():
    import networkx as nx
    from ccbench.support import build_support
    g = example()
    P = layout(g)
    sup = build_support(g)
    e = {tuple(sorted(x)) for x in map(tuple, g.edges())}
    two = [(u, v) for u, v in zip(sup.pu, sup.pv) if (min(u, v), max(u, v)) not in e]
    tot = g.n * (g.n - 1) // 2
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.0))
    draw_graph(axes[0], g, P, edge_color=INK2)
    axes[0].set_title(f"$G^+$: {g.n} vertices, {g.m} edges")
    ax = axes[1]
    for u, v in two:
        ax.plot(*P[[u, v]].T, color=ORANGE, lw=0.5, alpha=0.6, zorder=0)
    draw_graph(ax, g, P, edge_color=INK2)
    ax.set_title(f"support $P$: {g.m} edges + {len(two)} pairs at distance two")
    ax.plot([], [], color=INK2, lw=0.9, label="edge ($c_p=+1$)")
    ax.plot([], [], color=ORANGE, lw=1.2, label="distance two ($c_p=-1$)")
    fig.legend(loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.text(0.5, -0.09, f"The other {tot - g.m - len(two)} of the {tot} pairs are far "
             "and are fixed to separated.", ha="center", color=INK2, fontsize=7.5)
    save(fig, "fig_support")


# --------------------------------------------------------------------------
# 2. partition crossover
# --------------------------------------------------------------------------

def crossover_parts(g, a, b):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from ccbench.objective import normalize_labels
    a, b = normalize_labels(a), normalize_labels(b)
    ka, kb = a.max() + 1, b.max() + 1
    M = coo_matrix((np.ones(g.n), (a, b + ka)), shape=(ka + kb, ka + kb))
    _, comp = connected_components(M, directed=False)
    cx = comp[a]
    return a, b, cx


def comp_cost(g, lab, nodes):
    import ccbench as cc
    from ccbench.graph import induced_subgraph
    return cc.cost(induced_subgraph(g, nodes), lab[nodes])


def fig_crossover():
    import ccbench as cc
    from ccbench.memetic import partition_crossover
    from ccbench.pivot import pivot
    from ccbench.reduce import unit
    g = example()
    P = layout(g)
    a, b = pivot(g, 1), pivot(g, 5)
    child, _ = partition_crossover(unit(g), a, b)
    a, b, cx = crossover_parts(g, a, b)
    col = np.array([MUTED] * g.n, dtype=object)
    wins = {"a": [], "b": [], "tie": []}
    for c in np.unique(cx):
        nodes = np.flatnonzero(cx == c)
        ca, cb = comp_cost(g, a, nodes), comp_cost(g, b, nodes)
        key = "a" if ca < cb else "b" if cb < ca else "tie"
        wins[key].append((ca, cb))
        col[nodes] = {"a": BLUE, "b": ORANGE, "tie": MUTED}[key]
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 5.6))
    axes = axes.ravel()
    for ax, lab, title in ((axes[0], a, f"(a) parent $a$ (Pivot, seed 1): cost {cc.cost(g, a)}"),
                           (axes[1], b, f"(b) parent $b$ (Pivot, seed 5): cost {cc.cost(g, b)}")):
        hulls(ax, P, lab)
        draw_graph(ax, g, P, lab)
        ax.set_title(title)
    ax = axes[2]
    hulls(ax, P, cx, color="#e9e8e3", alpha=0.8)
    draw_graph(ax, g, P, None, node_color=list(col), edge_color=NEUTRAL)
    ax.set_title(f"(c) {len(np.unique(cx))} blocks of the overlap $a \\vee b$")
    ax = axes[3]
    hulls(ax, P, child)
    draw_graph(ax, g, P, child, node_color=list(col))
    ax.set_title(f"(d) child, better parent per block: cost {cc.cost(g, child)}")
    labs = [(BLUE, "block where $a$ is better (" + ", ".join(f"{x}<{y}" for x, y in wins["a"]) + ")"),
            (ORANGE, "block where $b$ is better (" + ", ".join(f"{y}<{x}" for x, y in wins["b"]) + ")"),
            (MUTED, f"tie ({len(wins['tie'])} blocks)")]
    for c, t in labs:
        axes[2].scatter([], [], s=34, c=c, edgecolors="white", label=t)
    axes[2].plot([], [], color=INK2, lw=1.1, label="edge inside a cluster")
    axes[2].plot([], [], color=NEUTRAL, lw=0.6, ls=(0, (2, 1.5)), label="cut edge")
    fig.legend(*axes[2].get_legend_handles_labels(), loc="lower center", ncol=2,
               frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    save(fig, "fig_crossover")


# --------------------------------------------------------------------------
# 3. gap map and gap-guided LNS
# --------------------------------------------------------------------------

def fig_gapmap():
    import ccbench as cc
    from ccbench.certiflip import block_bound
    from ccbench.lns import gap_lns, gap_map
    from ccbench.pivot import pivot
    g = example()
    P = layout(g)
    bd = block_bound(g, time_limit=20, block_size=200, seed=0)
    lb = bd.certified()
    lab0 = pivot(g, 7)
    lab1 = gap_lns(g, lab0, bd=bd, size=25, iters=60, time_limit=30, rng=0)
    vmax = None
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.1))
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    maps = []
    for lab in (lab0, lab1):
        vg, _, _ = gap_map(bd, lab)
        maps.append(vg)
    vmax = max(m.max() for m in maps)
    for ax, lab, vg, title in ((axes[0], lab0, maps[0], "Pivot"),
                               (axes[1], lab1, maps[1], "after gap-guided LNS")):
        hulls(ax, P, lab)
        e = g.edges()
        for u, v in e:
            inside = lab[u] == lab[v]
            ax.plot(*P[[u, v]].T, color=INK2 if inside else NEUTRAL,
                    lw=0.9 if inside else 0.5, zorder=1)
        sc = ax.scatter(P[:, 0], P[:, 1], s=30, c=vg, cmap=cmap, vmin=0, vmax=vmax,
                        edgecolors="white", linewidths=0.8, zorder=3)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_aspect("equal")
        c = cc.cost(g, lab)
        ax.set_title(f"{title}: cost {c}, gap {c - lb}" + (" (optimal)" if c == lb else ""))
    cb = fig.colorbar(sc, ax=axes, shrink=0.7, pad=0.02)
    cb.set_label("gap mass at the vertex", color=INK2)
    cb.outline.set_visible(False)
    fig.text(0.45, 0.08, f"certified lower bound {lb} (block dual ascent, 20 s)", ha="center",
             color=INK2, fontsize=7.5)
    save(fig, "fig_gapmap")


# --------------------------------------------------------------------------
# 4. convergence of PXMem against the final KaPoCE values (sequential runs)
# --------------------------------------------------------------------------

FIGS = {"support": fig_support, "crossover": fig_crossover, "gapmap": fig_gapmap}

if __name__ == "__main__":
    for k in (sys.argv[1:] or FIGS):
        FIGS[k]()
