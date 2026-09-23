"""Dual-gap guided large neighbourhood search with exact sub-MIPs.

Gap decomposition.  For multipliers y >= 0 with reduced costs r = c - A^T y
and any clustering with cut vector x^C (x_e = 1 iff the pair is separated),

    cost(C) - LB(y) = sum_e [ r_e x^C_e - min(0, r_e) ]  +  sum_i y_i (a_i x^C - b_i),

and every term is non-negative.  The pair terms vanish iff x^C agrees with the
sign of r, and the row terms iff every row with a positive multiplier is tight.
Distributing the terms to vertices gives a map of where C can still be
improved or where the bound is weak; a zero-gap region is certified.

Neighbourhoods.  A neighbourhood is a vertex set B grown around a high-gap
vertex.  With the clustering fixed outside B, the best re-clustering of B is a
weighted correlation clustering problem on B plus one super node per outside
cluster touching B, which we solve exactly as a MIP with lazily added
triangle inequalities.  A move is accepted iff the exact global cost drops.
"""
from __future__ import annotations

import time

import highspy
import numpy as np
import numba as nb

from .graph import Graph
from .objective import cost
from .localsearch import _compact


# --------------------------------------------------------------------------
# gap map
# --------------------------------------------------------------------------

def gap_map(bd, labels: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Per-vertex share of cost(C) - LB for a BlockDualBound ``bd``.

    Returns (vertex_gap, pair_part, row_part)."""
    sup = bd.sup
    lab = np.asarray(labels)
    xs = (lab[sup.pu] != lab[sup.pv]).astype(np.float64)
    r = bd.r
    ge = r * xs - np.minimum(r, 0.0)
    vg = np.zeros(bd.g.n)
    np.add.at(vg, sup.pu, 0.5 * ge)
    np.add.at(vg, sup.pv, 0.5 * ge)
    row_part = 0.0
    from .blockdual import _gather
    for ch in bd._chunks:
        sel = np.flatnonzero(ch["alive"])
        if len(sel) == 0:
            continue
        flat, rowof = _gather(ch["ptr"], sel)
        ax = np.zeros(len(sel))
        np.add.at(ax, rowof, ch["val"][flat] * xs[ch["idx"][flat]])
        slack = ch["y"][sel] * (ax - ch["b"][sel])
        slack = np.maximum(slack, 0.0)
        row_part += float(slack.sum())
        # distribute to the endpoints of the row's pairs
        lens = ch["ptr"][sel + 1] - ch["ptr"][sel]
        share = np.repeat(slack / np.maximum(1, 2 * lens), lens)
        np.add.at(vg, sup.pu[ch["idx"][flat]], share)
        np.add.at(vg, sup.pv[ch["idx"][flat]], share)
    return vg, float(ge.sum()), row_part


# --------------------------------------------------------------------------
# exact sub-problem
# --------------------------------------------------------------------------

@nb.njit(cache=True)
def _violated(k, pairidx, x, cap):
    """Violated triangle inequalities among k nodes; pairidx[i, j] = var id or -1 (x=1)."""
    out = np.empty((cap, 3), dtype=np.int64)
    typ = np.empty(cap, dtype=np.int8)
    c = 0
    for i in range(k):
        for j in range(i + 1, k):
            a = pairidx[i, j]
            xa = 1.0 if a < 0 else x[a]
            for l in range(j + 1, k):
                b = pairidx[i, l]
                cc = pairidx[j, l]
                xb = 1.0 if b < 0 else x[b]
                xc = 1.0 if cc < 0 else x[cc]
                # three inequalities: each side <= sum of other two
                if xa > xb + xc + 1e-6 or xb > xa + xc + 1e-6 or xc > xa + xb + 1e-6:
                    if c < cap:
                        out[c, 0] = a
                        out[c, 1] = b
                        out[c, 2] = cc
                        if xa > xb + xc + 1e-6:
                            typ[c] = 0
                        elif xb > xa + xc + 1e-6:
                            typ[c] = 1
                        else:
                            typ[c] = 2
                    c += 1
    m = min(c, cap)
    return out[:m], typ[:m], c


def _solve_sub(coef, const, pairidx, init_x=None, time_limit=30.0, max_rounds=50):
    k = pairidx.shape[0]
    N = len(coef)
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("mip_rel_gap", 0.0)
    h.addVars(N, np.zeros(N), np.ones(N))
    h.changeColsCost(N, np.arange(N, dtype=np.int32), coef)
    h.changeColsIntegrality(N, np.arange(N, dtype=np.int32),
                            np.array([highspy.HighsVarType.kInteger] * N))
    t0 = time.time()
    x = None
    for _ in range(max_rounds):
        h.setOptionValue("time_limit", max(1.0, time_limit - (time.time() - t0)))
        if init_x is not None:
            try:
                sol = highspy.HighsSolution()
                sol.col_value = list(init_x)
                h.setSolution(sol)
            except Exception:
                pass
        h.run()
        if h.getModelStatus() not in (highspy.HighsModelStatus.kOptimal,
                                      highspy.HighsModelStatus.kTimeLimit,
                                      highspy.HighsModelStatus.kInterrupt):
            return None
        x = np.round(np.asarray(h.getSolution().col_value))
        tri, typ, found = _violated(k, pairidx, x, 20000)
        if found == 0:
            return x
        if time.time() - t0 > time_limit:
            return None
        # add rows: x_big <= x_1 + x_2 ; fixed pairs (id -1, x = 1) move to rhs
        starts, index, value, lower = [], [], [], []
        nnz = 0
        for (a, b, c), t in zip(tri, typ):
            trip = [a, b, c]
            big = trip[t]
            rest = [trip[j] for j in range(3) if j != t]
            # rest_1 + rest_2 - big >= 0
            idxs, vals, lo = [], [], 0.0
            for q in rest:
                if q < 0:
                    lo -= 1.0
                else:
                    idxs.append(q)
                    vals.append(1.0)
            if big < 0:
                lo += 1.0
            else:
                idxs.append(big)
                vals.append(-1.0)
            if not idxs:
                continue
            starts.append(nnz)
            index += idxs
            value += vals
            lower.append(lo)
            nnz += len(idxs)
        if not starts:
            return x
        h.addRows(len(lower), np.array(lower), np.full(len(lower), highspy.kHighsInf), nnz,
                  np.array(starts, dtype=np.int32), np.array(index, dtype=np.int32),
                  np.array(value))
    return None


@nb.njit(cache=True)
def _within2_nb(indptr, indices, u, v):
    """True iff u and v are adjacent or share a neighbour in G+."""
    i = indptr[u]
    j = indptr[v]
    while i < indptr[u + 1] and j < indptr[v + 1]:
        a = indices[i]
        b = indices[j]
        if a == b or a == v or b == u:
            return True
        if a < b:
            i += 1
        else:
            j += 1
    while i < indptr[u + 1]:
        if indices[i] == v:
            return True
        i += 1
    while j < indptr[v + 1]:
        if indices[j] == u:
            return True
        j += 1
    return False


def _within2(g: Graph, u: int, v: int) -> bool:
    return bool(_within2_nb(g.indptr, g.indices, u, v))


def reoptimize_block(g: Graph, labels: np.ndarray, B: np.ndarray, time_limit: float = 20.0,
                     max_join: int = 8, max_super: int = 500):
    """Best re-clustering of the vertices B with everything else fixed, among
    re-clusterings in which no pair at G+-distance >= 3 shares a cluster.

    Returns new labels (or None if the sub-MIP failed)."""
    lab = np.asarray(labels).copy()
    B = np.asarray(B, dtype=np.int64)
    inB = np.zeros(g.n, dtype=bool)
    inB[B] = True
    posB = np.full(g.n, -1, dtype=np.int64)
    posB[B] = np.arange(len(B))
    # outside clusters touching B via positive edges
    csize = np.bincount(lab, minlength=lab.max() + 1)
    inside_count = np.bincount(lab[B], minlength=len(csize))
    outer_size = csize - inside_count
    wK = {}
    adjB = []
    for i, v in enumerate(B):
        nb_ = g.indices[g.indptr[v]:g.indptr[v + 1]]
        adjB.append(nb_)
        for u in nb_:
            if not inB[u]:
                key = (i, int(lab[u]))
                wK[key] = wK.get(key, 0) + 1
    # restrict join options: each vertex may join its (at most) ``max_join`` best
    # outside clusters, and at most ``max_super`` super nodes overall; other
    # options are fixed apart (their cost is a constant of the sub-problem).
    n_super = len({K for (_, K) in wK})
    per_v = {}
    for (i, K), w in (wK.items() if n_super > max_super else ()):
        per_v.setdefault(i, []).append((w, K))
    keep = {} if n_super > max_super else dict(wK)
    for i, lst in per_v.items():
        lst.sort(reverse=True)
        for w, K in lst[:max_join]:
            keep[(i, K)] = w
    tot = {}
    for (i, K), w in keep.items():
        tot[K] = tot.get(K, 0) + w
    if len(tot) > max_super:
        allowed = set(sorted(tot, key=lambda K: -tot[K])[:max_super])
        keep = {key: w for key, w in keep.items() if key[1] in allowed}
    wK = keep  # dropped options only add a constant to the sub-problem
    Ks = sorted({K for (_, K) in wK})
    kpos = {K: len(B) + j for j, K in enumerate(Ks)}
    k = len(B) + len(Ks)
    pairidx = -np.ones((k, k), dtype=np.int64)
    coef, cst = [], 0.0
    E = {}
    for i, v in enumerate(B):
        for u in adjB[i]:
            j = posB[u]
            if j > i:
                E[(i, j)] = True
    nv = 0
    for i in range(len(B)):
        for j in range(i + 1, len(B)):
            if (i, j) not in E and not _within2(g, int(B[i]), int(B[j])):
                continue  # far pair: must be separated
            pairidx[i, j] = pairidx[j, i] = nv
            nv += 1
            if (i, j) in E:
                coef.append(1.0)           # cut positive
            else:
                coef.append(-1.0)          # together costs 1 - x
                cst += 1.0
    # a vertex may join an outside cluster only if it is within G+-distance 2 of
    # all of its members (P-feasibility: optimal clusters have diameter <= 2)
    # members of the touched outside clusters (one pass)
    if Ks:
        outside = np.flatnonzero(np.isin(lab, Ks) & ~inB)
        members = {K: [] for K in Ks}
        for u in outside:
            members[int(lab[u])].append(int(u))
    for (i, K), w in wK.items():
        v = int(B[i])
        if any(not _within2_nb(g.indptr, g.indices, v, u) for u in members[K]):
            continue
        j = kpos[K]
        s = float(outer_size[K])
        pairidx[i, j] = pairidx[j, i] = nv
        nv += 1
        # w x + (s - w)(1 - x)
        coef.append(2.0 * w - s)
        cst += s - w
    coef = np.array(coef)
    # incumbent as a start
    x0 = np.ones(nv)
    for i in range(len(B)):
        for j in range(i + 1, len(B)):
            if pairidx[i, j] >= 0:
                x0[pairidx[i, j]] = float(lab[B[i]] != lab[B[j]])
    for (i, K) in wK:
        if pairidx[i, kpos[K]] >= 0:
            x0[pairidx[i, kpos[K]]] = float(lab[B[i]] != K)
    x = _solve_sub(coef, cst, pairidx, x0, time_limit)
    if x is None:
        return None
    # decode: union of x = 0 pairs
    parent = np.arange(k)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for i in range(k):
        for j in range(i + 1, k):
            q = pairidx[i, j]
            if q >= 0 and x[q] < 0.5:
                parent[find(i)] = find(j)
    root_label = {}
    for K in Ks:
        root_label[find(kpos[K])] = K
    nxt = lab.max() + 1
    for i, v in enumerate(B):
        r0 = find(i)
        if r0 not in root_label:
            root_label[r0] = nxt
            nxt += 1
        lab[v] = root_label[r0]
    return lab


def grow_block(g: Graph, seed: int, size: int, score: np.ndarray, rng) -> np.ndarray:
    """BFS-like growth preferring high-score neighbours."""
    chosen = [seed]
    seen = {seed}
    frontier = {}
    for u in g.indices[g.indptr[seed]:g.indptr[seed + 1]]:
        frontier[int(u)] = score[u]
    while len(chosen) < size and frontier:
        u = max(frontier, key=lambda z: frontier[z] + 1e-6 * rng.random())
        del frontier[u]
        chosen.append(u)
        seen.add(u)
        for w in g.indices[g.indptr[u]:g.indptr[u + 1]]:
            w = int(w)
            if w not in seen and w not in frontier:
                frontier[w] = score[w]
    return np.array(chosen, dtype=np.int64)


def touching_gap(bd, labels: np.ndarray, B: np.ndarray) -> float:
    """Gap mass of all terms (pairs and rows) touching the vertex set B.

    Local certificate: no re-clustering of B (outside fixed) can decrease the
    cost by more than this amount, because every term not touching B is
    unchanged and every term is non-negative."""
    sup = bd.sup
    lab = np.asarray(labels)
    inB = np.zeros(bd.g.n if bd.g is not None else len(lab), dtype=bool)
    inB[B] = True
    touch = inB[sup.pu] | inB[sup.pv]
    xs = (lab[sup.pu[touch]] != lab[sup.pv[touch]]).astype(np.float64)
    r = bd.r[touch]
    tot = float((r * xs - np.minimum(r, 0.0)).sum())
    from .blockdual import _gather
    xfull = None
    for ch in bd._chunks:
        sel = np.flatnonzero(ch["alive"])
        if len(sel) == 0:
            continue
        flat, rowof = _gather(ch["ptr"], sel)
        pids = ch["idx"][flat]
        t_row = np.zeros(len(sel), dtype=bool)
        np.logical_or.at(t_row, rowof, inB[sup.pu[pids]] | inB[sup.pv[pids]])
        if not t_row.any():
            continue
        if xfull is None:
            xfull = (lab[sup.pu] != lab[sup.pv]).astype(np.float64)
        ax = np.zeros(len(sel))
        np.add.at(ax, rowof, ch["val"][flat] * xfull[pids])
        slack = np.maximum(ch["y"][sel] * (ax - ch["b"][sel]), 0.0)
        tot += float(slack[t_row].sum())
    return tot


def cluster_block(g: Graph, labels: np.ndarray, seed: int, size: int, score: np.ndarray):
    """Union of the seed's cluster and neighbouring clusters (highest total score
    first), truncated to ``size`` vertices by score."""
    lab = np.asarray(labels)
    members = {}
    order = np.argsort(lab, kind="stable")
    # cheap member lookup for the touched clusters only
    cand = {int(lab[seed])}
    for u in g.indices[g.indptr[seed]:g.indptr[seed + 1]]:
        cand.add(int(lab[u]))
    mask = np.isin(lab, list(cand))
    verts = np.flatnonzero(mask)
    for v in verts:
        members.setdefault(int(lab[v]), []).append(int(v))
    ranked = sorted(members, key=lambda c: -sum(score[members[c]]))
    ranked.remove(int(lab[seed]))
    ranked = [int(lab[seed])] + ranked
    out = []
    for c in ranked:
        out += members[c]
        if len(out) >= size:
            break
    out = np.array(out, dtype=np.int64)
    if len(out) > size:
        keep = np.argsort(-score[out])[:size]
        out = out[keep]
        if seed not in set(out.tolist()):
            out[-1] = seed
    return out


def gap_lns(g: Graph, labels: np.ndarray, bd=None, size: int = 60, iters: int = 200,
            time_limit: float = 300.0, sub_time: float = 10.0, rng=None, verbose=False,
            modes=("bfs", "cluster"), stats: dict | None = None):
    """Dual-gap guided LNS.  Without a bound ``bd``, seeds are chosen by the
    local disagreement count.  Neighbourhoods whose touching gap is < 1 are
    certified locally optimal and skipped."""
    rng = np.random.default_rng(rng)
    lab = _compact(np.asarray(labels))[0].copy()
    cur = cost(g, lab)
    t0 = time.time()
    tabu = np.zeros(g.n, dtype=np.int64)
    it = 0
    st = {"tried": 0, "improved": 0, "skipped": 0, "failed": 0}
    score = None
    while it < iters and time.time() - t0 < time_limit:
        if score is None or st["improved"] or it % 10 == 0:
            if bd is not None:
                score, _, _ = gap_map(bd, lab)
            else:
                src = np.repeat(np.arange(g.n), g.degrees)
                cut = lab[src] != lab[g.indices]
                score = np.bincount(src[cut], minlength=g.n).astype(float)
        eff = score - 1e3 * (tabu > it)
        if eff.max() <= 1e-9:
            break
        order = np.argsort(-eff)
        seed = int(order[rng.integers(0, min(20, g.n))])
        mode = modes[it % len(modes)]
        if mode == "bfs":
            B = grow_block(g, seed, size, score, rng)
        else:
            B = cluster_block(g, lab, seed, size, score)
        tabu[B] = it + 10
        it += 1
        if bd is not None and touching_gap(bd, lab, B) < 1.0 - 1e-6:
            st["skipped"] += 1
            continue
        st["tried"] += 1
        new = reoptimize_block(g, lab, B, sub_time)
        if new is None:
            st["failed"] += 1
            continue
        c = cost(g, new)
        if c < cur:
            lab = _compact(new)[0].copy()
            st["improved"] += 1
            if verbose:
                print(f"  lns it {it} ({mode}): {cur} -> {c} ({time.time() - t0:.1f}s)", flush=True)
            cur = c
    if stats is not None:
        stats.update(st)
    return lab
