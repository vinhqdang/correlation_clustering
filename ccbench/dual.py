"""Lower bounds for correlation clustering on complete signed graphs.

Lagrangian (dual-decomposition) bound
-------------------------------------
Write the problem over the distance-2 support P with cut indicators y_e
(y_e = 1: pair separated).  Far pairs are cut in every optimal solution, so

    OPT = |N2| + min  sum_e c_e y_e ,  c_e = +1 (e in E+), -1 (e in N2),

over cut vectors of clusterings.  For any family T of triangles in P, each
triangle's restriction lies in F = {000, 110, 101, 011, 111}.  For arbitrary
multipliers lambda_{t,e} with theta_e = c_e - sum_{t ni e} lambda_{t,e},

    LB(lambda) = |N2| + sum_e min(0, theta_e) + sum_t min_{y in F} <lambda_t, y>

is a valid lower bound on OPT (weak duality; it equals the value of the
standard LP restricted to the triangle inequalities of T at the optimum over
lambda).  We maximise LB by block-coordinate ascent: every triangle update is
monotone, so the reported bound never decreases and is valid at any time.

T contains every triangle of P with at least two positive pairs, i.e. all
bad triangles (open wedges) and all positive triangles of G+.

Packing bound
-------------
A set of pair-disjoint bad triangles certifies one disagreement each.
"""
from __future__ import annotations

import time

import numpy as np
import numba as nb

from .graph import Graph
from .support import Support, build_support, wedge_triangles


@nb.njit(cache=True, inline="always")
def _tmin(x, y, z):
    m = 0.0
    if x + y < m:
        m = x + y
    if x + z < m:
        m = x + z
    if y + z < m:
        m = y + z
    if x + y + z < m:
        m = x + y + z
    return m


@nb.njit(cache=True, inline="always")
def _mm(x, y, z):
    """Min-marginal of the first coordinate: min_{y1=1} - min_{y1=0}."""
    m1 = y
    if z < m1:
        m1 = z
    if y + z < m1:
        m1 = y + z
    m0 = 0.0
    if y + z < m0:
        m0 = y + z
    return x + m1 - m0


@nb.njit(cache=True)
def _sweep(ta, tb, tc, la, lb, lc, theta, cnt, order):
    for k in range(order.shape[0]):
        t = order[k]
        a = ta[t]
        b = tb[t]
        c = tc[t]
        # pull a share of the unary costs into the triangle
        da = theta[a] / cnt[a]
        db = theta[b] / cnt[b]
        dc = theta[c] / cnt[c]
        x = la[t] + da
        y = lb[t] + db
        z = lc[t] + dc
        theta[a] -= da
        theta[b] -= db
        theta[c] -= dc
        # push back min-marginals (sequential, each step keeps the bound)
        m = _mm(x, y, z)
        x -= m
        theta[a] += m
        m = _mm(y, x, z)
        y -= m
        theta[b] += m
        m = _mm(z, x, y)
        z -= m
        theta[c] += m
        la[t] = x
        lb[t] = y
        lc[t] = z


@nb.njit(cache=True)
def _bound(ta, tb, tc, la, lb, lc, theta, const):
    s = const
    for e in range(theta.shape[0]):
        if theta[e] < 0:
            s += theta[e]
    for t in range(ta.shape[0]):
        s += _tmin(la[t], lb[t], lc[t])
    return s


@nb.njit(cache=True)
def _counts(ta, tb, tc, npairs):
    cnt = np.ones(npairs, dtype=np.float64)
    for t in range(ta.shape[0]):
        cnt[ta[t]] += 1.0
        cnt[tb[t]] += 1.0
        cnt[tc[t]] += 1.0
    return cnt


class LagrangianBound:
    """Block-coordinate ascent on the triangle dual decomposition."""

    def __init__(self, g: Graph, sup: Support | None = None, tris=None):
        self.g = g
        self.sup = sup if sup is not None else build_support(g)
        if tris is None:
            tris = wedge_triangles(self.sup)
        idt = np.int32 if self.sup.npairs < 2**31 - 1 else np.int64
        self.ta, self.tb, self.tc = (np.ascontiguousarray(x, dtype=idt) for x in tris)
        T = self.ta.shape[0]
        self.la = np.zeros(T)
        self.lb = np.zeros(T)
        self.lc = np.zeros(T)
        self.theta = np.empty(self.sup.npairs)
        self.theta[: self.sup.npos] = 1.0
        self.theta[self.sup.npos:] = -1.0
        self.const = float(self.sup.nneg)
        self.cnt = _counts(self.ta, self.tb, self.tc, self.sup.npairs)
        self.history = []

    @property
    def ntriangles(self) -> int:
        return int(self.ta.shape[0])

    def bound(self) -> float:
        return _bound(self.ta, self.tb, self.tc, self.la, self.lb, self.lc,
                      self.theta, self.const)

    def run(self, sweeps: int = 100, tol: float = 1e-4, time_limit: float = 1e9,
            rng=None, check_every: int = 5) -> float:
        rng = np.random.default_rng(rng)
        order = np.arange(self.ntriangles, dtype=np.int64)
        t0 = time.time()
        last = self.bound()
        self.history.append((0, 0.0, last))
        for it in range(1, sweeps + 1):
            _sweep(self.ta, self.tb, self.tc, self.la, self.lb, self.lc,
                   self.theta, self.cnt, order)
            if it % check_every == 0 or it == sweeps:
                cur = self.bound()
                self.history.append((it, time.time() - t0, cur))
                if cur - last < tol * max(1.0, abs(cur)):
                    break
                last = cur
            if time.time() - t0 > time_limit:
                break
        return self.bound()

    def certified_bound(self) -> int:
        """Integer lower bound (OPT is integral)."""
        return int(np.ceil(self.bound() - 1e-6))

    def reduced_costs(self) -> np.ndarray:
        """theta plus the triangle min-marginals: a per-pair 'cut preference'.

        Negative values mean the dual prefers separating the pair."""
        return _reduced(self.ta, self.tb, self.tc, self.la, self.lb, self.lc,
                        self.theta.copy())


@nb.njit(cache=True)
def _reduced(ta, tb, tc, la, lb, lc, r):
    for t in range(ta.shape[0]):
        x = la[t]
        y = lb[t]
        z = lc[t]
        r[ta[t]] += _mm(x, y, z)
        r[tb[t]] += _mm(y, x, z)
        r[tc[t]] += _mm(z, x, y)
    return r


@nb.njit(cache=True)
def _greedy_packing(n, indptr, indices, eid, n2ptr, n2idx, n2id, npairs, order):
    used = np.zeros(npairs, dtype=np.bool_)
    cnt = 0
    for k in range(n):
        w = order[k]
        s = indptr[w]
        e = indptr[w + 1]
        for p in range(s, e):
            if used[eid[p]]:
                continue
            u = indices[p]
            for q in range(p + 1, e):
                if used[eid[q]]:
                    continue
                v = indices[q]
                # negative pair uv?
                lo = indptr[u]
                hi = indptr[u + 1]
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if indices[mid] < v:
                        lo = mid + 1
                    else:
                        hi = mid
                if lo < indptr[u + 1] and indices[lo] == v:
                    continue
                lo = n2ptr[u]
                hi = n2ptr[u + 1]
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if n2idx[mid] < v:
                        lo = mid + 1
                    else:
                        hi = mid
                r = n2id[lo]
                if used[r]:
                    continue
                used[eid[p]] = True
                used[eid[q]] = True
                used[r] = True
                cnt += 1
                break
    return cnt


def packing_bound(g: Graph, sup: Support | None = None, rng=None) -> int:
    """Maximal pair-disjoint packing of bad triangles (low-degree centres first)."""
    sup = sup if sup is not None else build_support(g)
    order = np.argsort(g.degrees, kind="stable").astype(np.int64)
    return int(_greedy_packing(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx,
                               sup.n2id, sup.npairs, order))


# --------------------------------------------------------------------------
# (1 - eps)-approximate fractional packing of bad triangles
# --------------------------------------------------------------------------
#
# The covering LP  min sum_{e in E+} x_e + sum_{e in N2} z_e  s.t.
# x_a + x_b + z_c >= 1 for every bad triangle (a, b positive, c negative)
# is a relaxation of CC; its dual is the fractional packing
#     max sum_t y_t   s.t.  sum_{t ni e} y_t <= 1,  y >= 0.
# We solve it with the Garg-Koenemann / Fleischer multiplicative-weights
# scheme.  The returned bound is sum_t y_t divided by the maximum load, so it
# is a valid lower bound whatever the parameters.  The final lengths give a
# feasible covering solution, i.e. an upper bound on the covering LP value.

@nb.njit(cache=True)
def _mwu_packing(ta, tb, tc, npairs, eps, delta, max_phases):
    T = ta.shape[0]
    ln = np.full(npairs, delta)
    y = np.zeros(T)
    active = np.arange(T)
    na = T
    alpha = 3.0 * delta
    f = 1.0 + eps
    phases = 0
    while alpha < 1.0 and phases < max_phases and na > 0:
        thr = alpha * f
        if thr > 1.0:
            thr = 1.0
        k = 0
        for i in range(na):
            t = active[i]
            a = ta[t]
            b = tb[t]
            c = tc[t]
            L = ln[a] + ln[b] + ln[c]
            while L < thr:
                y[t] += 1.0
                ln[a] *= f
                ln[b] *= f
                ln[c] *= f
                L = ln[a] + ln[b] + ln[c]
            if L < 1.0:
                active[k] = t
                k += 1
        na = k
        alpha = thr
        phases += 1
    load = np.zeros(npairs)
    for t in range(T):
        load[ta[t]] += y[t]
        load[tb[t]] += y[t]
        load[tc[t]] += y[t]
    mx = 0.0
    for e in range(npairs):
        if load[e] > mx:
            mx = load[e]
    tot = 0.0
    for t in range(T):
        tot += y[t]
    # covering solution from lengths
    mn = 1e300
    for t in range(T):
        L = ln[ta[t]] + ln[tb[t]] + ln[tc[t]]
        if L < mn:
            mn = L
    return y, mx, tot, ln, mn, phases


class PackingBound:
    """(1-eps)-approximate fractional bad-triangle packing (covering-LP dual)."""

    def __init__(self, g: Graph, sup: Support | None = None, tris=None):
        self.g = g
        self.sup = sup if sup is not None else build_support(g)
        if tris is None:
            tris = wedge_triangles(self.sup)
        ta, tb, tc = tris
        bad = tc >= self.sup.npos
        idt = np.int32 if self.sup.npairs < 2**31 - 1 else np.int64
        self.ta = np.ascontiguousarray(ta[bad], dtype=idt)
        self.tb = np.ascontiguousarray(tb[bad], dtype=idt)
        self.tc = np.ascontiguousarray(tc[bad], dtype=idt)

    def run(self, eps: float = 0.1, delta: float | None = None, max_phases: int = 100000):
        M = self.sup.npairs
        if delta is None:
            delta = (1 + eps) * ((1 + eps) * max(M, 2)) ** (-1.0 / eps)
        t0 = time.time()
        y, mx, tot, ln, mn, ph = _mwu_packing(self.ta, self.tb, self.tc, M, eps, delta,
                                              max_phases)
        self.time = time.time() - t0
        self.phases = ph
        self.y = y / max(mx, 1e-300)
        self.lower = tot / max(mx, 1e-300) if tot > 0 else 0.0
        # covering upper bound (only pairs that lie in some bad triangle matter)
        used = np.zeros(M, dtype=bool)
        used[self.ta] = used[self.tb] = used[self.tc] = True
        self.cover = np.where(used, np.minimum(1.0, ln / mn), 0.0) if mn > 0 else None
        self.upper = float(self.cover.sum()) if self.cover is not None else float("inf")
        return self.lower


# --------------------------------------------------------------------------
# fractional packing of induced stars
# --------------------------------------------------------------------------
#
# An induced star (v; T) with T a set of k >= 2 pairwise non-adjacent positive
# neighbours of v forces at least k - 1 disagreements among its k + k(k-1)/2
# pairs (the star inequality in "mistake" form).  A fractional packing
#     max sum_S (k_S - 1) y_S   s.t.   sum_{S ni e} y_S <= 1
# is therefore a lower bound on OPT; bad triangles are the stars with k = 2.
# We run the Garg-Koenemann scheme with a greedy column oracle per centre and
# normalise by the maximum load at the end, so the bound is always valid.

@nb.njit(cache=True)
def _best_star(v, indptr, indices, eid, ptr, idx, pid, ln, cap_k, cand, candl, T, Tp):
    """Greedy star of minimum length/value ratio centred at v.

    Returns (k, total_length); T[:k] leaves, Tp[:k + k(k-1)/2] pair ids."""
    d = indptr[v + 1] - indptr[v]
    if d < 2:
        return 0, 0.0
    for i in range(d):
        p = indptr[v] + i
        cand[i] = indices[p]
        candl[i] = ln[eid[p]]
    order = np.argsort(candl[:d])
    k = 0
    npair = 0
    L = 0.0
    for oi in range(d):
        i = order[oi]
        t = cand[i]
        # independence and added length
        ok = True
        add = candl[i]
        for j in range(k):
            lo = indptr[t]
            hi = indptr[t + 1]
            while lo < hi:
                mid = (lo + hi) >> 1
                if indices[mid] < T[j]:
                    lo = mid + 1
                else:
                    hi = mid
            if lo < indptr[t + 1] and indices[lo] == T[j]:
                ok = False
                break
            lo = ptr[t]
            hi = ptr[t + 1]
            while lo < hi:
                mid = (lo + hi) >> 1
                if idx[mid] < T[j]:
                    lo = mid + 1
                else:
                    hi = mid
            add += ln[pid[lo]]
        if not ok:
            continue
        if k >= 2 and (L + add) / k >= L / (k - 1):
            continue
        # accept t: record pair ids (vt and tt' for t' in T)
        Tp[npair] = eid[indptr[v] + i]
        npair += 1
        for j in range(k):
            lo = ptr[t]
            hi = ptr[t + 1]
            while lo < hi:
                mid = (lo + hi) >> 1
                if idx[mid] < T[j]:
                    lo = mid + 1
                else:
                    hi = mid
            Tp[npair] = pid[lo]
            npair += 1
        T[k] = t
        k += 1
        L += add
        if k == cap_k:
            break
    if k < 2:
        return 0, 0.0
    return k, L


@nb.njit(cache=True)
def _star_mwu(n, indptr, indices, eid, ptr, idx, pid, npairs, eps, delta, max_phases, cap_k,
              max_aug):
    ln = np.full(npairs, delta)
    load = np.zeros(npairs)
    maxd = 0
    for v in range(n):
        if indptr[v + 1] - indptr[v] > maxd:
            maxd = indptr[v + 1] - indptr[v]
    cand = np.empty(maxd + 1, dtype=np.int64)
    candl = np.empty(maxd + 1, dtype=np.float64)
    T = np.empty(cap_k + 1, dtype=np.int64)
    Tp = np.empty(cap_k + cap_k * (cap_k + 1) // 2 + 2, dtype=np.int64)
    f = 1.0 + eps
    alpha = 3.0 * delta
    value = 0.0
    phases = 0
    active = np.ones(n, dtype=np.bool_)
    while alpha < 1.0 and phases < max_phases:
        thr = min(1.0, alpha * f)
        any_active = False
        for v in range(n):
            if not active[v]:
                continue
            aug = 0
            while aug < max_aug:
                k, L = _best_star(v, indptr, indices, eid, ptr, idx, pid, ln, cap_k, cand,
                                  candl, T, Tp)
                if k < 2:
                    active[v] = False
                    break
                r = L / (k - 1)
                if r >= 1.0:
                    active[v] = False
                    break
                if r >= thr:
                    break
                npair = k + k * (k - 1) // 2
                for j in range(npair):
                    e = Tp[j]
                    load[e] += 1.0
                    ln[e] *= f
                value += k - 1
                aug += 1
            if active[v]:
                any_active = True
        if not any_active:
            break
        alpha = thr
        phases += 1
    mx = 0.0
    for e in range(npairs):
        if load[e] > mx:
            mx = load[e]
    return value, mx, phases


class StarPackingBound:
    """Fractional induced-star packing lower bound (generalises bad triangles)."""

    def __init__(self, g: Graph, sup: Support | None = None, pgraph=None):
        from .lp import _pgraph
        self.g = g
        self.sup = sup if sup is not None else build_support(g)
        s = self.sup
        self.pgraph = pgraph if pgraph is not None else _pgraph(
            g.n, g.indptr, g.indices, s.eid, s.n2ptr, s.n2idx, s.n2id)

    def run(self, eps: float = 0.1, delta: float = 1e-4, max_phases: int = 100000,
            cap_k: int = 64, max_aug: int = 1 << 30):
        g, s = self.g, self.sup
        ptr, idx, pid = self.pgraph
        t0 = time.time()
        value, mx, ph = _star_mwu(g.n, g.indptr, g.indices, s.eid, ptr, idx, pid, s.npairs,
                                  eps, delta, max_phases, cap_k, max_aug)
        self.time = time.time() - t0
        self.phases = ph
        self.lower = value / mx if mx > 0 else 0.0
        return self.lower


@nb.njit(cache=True)
def _greedy_packing_mask(n, indptr, indices, eid, n2ptr, n2idx, n2id, npairs, order):
    """Like _greedy_packing but returns the mask of pairs used by the packing."""
    used = np.zeros(npairs, dtype=np.bool_)
    for k in range(n):
        w = order[k]
        s = indptr[w]
        e = indptr[w + 1]
        for p in range(s, e):
            if used[eid[p]]:
                continue
            u = indices[p]
            for q in range(p + 1, e):
                if used[eid[q]] or used[eid[p]]:
                    continue
                v = indices[q]
                lo = indptr[u]
                hi = indptr[u + 1]
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if indices[mid] < v:
                        lo = mid + 1
                    else:
                        hi = mid
                if lo < indptr[u + 1] and indices[lo] == v:
                    continue
                lo = n2ptr[u]
                hi = n2ptr[u + 1]
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if n2idx[mid] < v:
                        lo = mid + 1
                    else:
                        hi = mid
                r = n2id[lo]
                if used[r]:
                    continue
                used[eid[p]] = True
                used[eid[q]] = True
                used[r] = True
    return used


def match_flip_pivot(g: Graph, sup: Support | None = None, rng=None):
    """MatchFlipPivot (Veldt, ICML 2022, Algorithm 3): maximal pair-disjoint set of
    open wedges (a maximal matching of the open-wedge hypergraph), flip every
    pair of the matched wedges, and run Pivot on the flipped graph.  Returns
    (labels, lower bound = number of matched wedges)."""
    from .graph import from_edges
    from .pivot import pivot
    rng = np.random.default_rng(rng)
    sup = sup if sup is not None else build_support(g)
    order = rng.permutation(g.n).astype(np.int64)
    used = _greedy_packing_mask(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx,
                                sup.n2id, sup.npairs, order)
    pos = np.arange(sup.npairs) < sup.npos
    keep = pos ^ used
    h = from_edges(g.n, np.stack([sup.pu[keep], sup.pv[keep]], axis=1))
    return pivot(h, rng), int(used.sum() // 3)


# --------------------------------------------------------------------------
# greedy integral packing of induced stars (sparse)
# --------------------------------------------------------------------------

@nb.njit(cache=True)
def _star_pack(n, indptr, indices, eid, ptr, idx, pid, npairs, order, leaf_rank):
    """Pair-disjoint packing of induced stars (centre v, >= 2 pairwise
    non-adjacent leaves).  A star with k leaves certifies k - 1 mistakes.

    Returns (value, row_ptr, row_pairs, row_sizes): each star's pair ids
    (centre pairs first, then leaf pairs) and its number of leaves."""
    used = np.zeros(npairs, dtype=np.bool_)
    cap = npairs + 1
    rp = np.zeros(npairs // 3 + 2, dtype=np.int64)
    rpairs = np.empty(cap, dtype=np.int64)
    rk = np.empty(npairs // 3 + 1, dtype=np.int64)
    ns = 0
    nz = 0
    maxd = 0
    for v in range(n):
        if indptr[v + 1] - indptr[v] > maxd:
            maxd = indptr[v + 1] - indptr[v]
    T = np.empty(maxd + 1, dtype=np.int64)
    Tc = np.empty(maxd + 1, dtype=np.int64)
    inner = np.empty(maxd * 4 + 16, dtype=np.int64)
    value = 0
    for oi in range(n):
        v = order[oi]
        d = indptr[v + 1] - indptr[v]
        if d < 2:
            continue
        # scan leaves in the order given by leaf_rank (e.g. low degree first)
        cand = np.empty(d, dtype=np.int64)
        for i in range(d):
            cand[i] = indptr[v] + i
        keys = np.empty(d, dtype=np.int64)
        for i in range(d):
            keys[i] = leaf_rank[indices[cand[i]]]
        o = np.argsort(keys)
        while True:
            k = 0
            ni = 0
            for oj in range(d):
                p = cand[o[oj]]
                if used[eid[p]]:
                    continue
                t = indices[p]
                ok = True
                start_ni = ni
                for j in range(k):
                    s = T[j]
                    # t and s must be non-adjacent
                    lo = indptr[t]
                    hi = indptr[t + 1]
                    while lo < hi:
                        mid = (lo + hi) >> 1
                        if indices[mid] < s:
                            lo = mid + 1
                        else:
                            hi = mid
                    if lo < indptr[t + 1] and indices[lo] == s:
                        ok = False
                        break
                    lo = ptr[t]
                    hi = ptr[t + 1]
                    while lo < hi:
                        mid = (lo + hi) >> 1
                        if idx[mid] < s:
                            lo = mid + 1
                        else:
                            hi = mid
                    e2 = pid[lo]
                    if used[e2]:
                        ok = False
                        break
                    if ni >= inner.shape[0]:
                        ok = False
                        break
                    inner[ni] = e2
                    ni += 1
                if not ok:
                    ni = start_ni
                    continue
                T[k] = t
                Tc[k] = eid[p]
                k += 1
            if k < 2:
                break
            for j in range(k):
                used[Tc[j]] = True
                rpairs[nz] = Tc[j]
                nz += 1
            for j in range(ni):
                used[inner[j]] = True
                rpairs[nz] = inner[j]
                nz += 1
            rk[ns] = k
            ns += 1
            rp[ns] = nz
            value += k - 1
    return value, rp[:ns + 1], rpairs[:nz], rk[:ns]


def star_packing(g: Graph, sup: Support | None = None, pgraph=None, order: str = "degree",
                 rng=None):
    """Greedy integral star packing.  Returns (value, rows) with rows usable as
    initial multipliers (y = 1) in the block dual."""
    from .lp import _pgraph
    sup = sup if sup is not None else build_support(g)
    if pgraph is None:
        pgraph = _pgraph(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id)
    ptr, idx, pid = pgraph
    deg = g.degrees
    rng = np.random.default_rng(rng)
    if order == "degree":
        ordv = np.lexsort((rng.random(g.n), deg)).astype(np.int64)
    elif order == "degree-desc":
        ordv = np.lexsort((rng.random(g.n), -deg)).astype(np.int64)
    else:
        ordv = rng.permutation(g.n).astype(np.int64)
    leaf_rank = np.argsort(np.argsort(deg + rng.random(g.n))).astype(np.int64)
    return _star_pack(g.n, g.indptr, g.indices, sup.eid, ptr, idx, pid, sup.npairs, ordv,
                      leaf_rank)


# --------------------------------------------------------------------------
# ruin-and-recreate local search for the star packing
# --------------------------------------------------------------------------

@nb.njit(cache=True)
def _find_pair(indptr, indices, ptr, idx, pid, t, s):
    """(adjacent?, pair id of ts in P or -1)."""
    lo = indptr[t]
    hi = indptr[t + 1]
    while lo < hi:
        mid = (lo + hi) >> 1
        if indices[mid] < s:
            lo = mid + 1
        else:
            hi = mid
    adj = lo < indptr[t + 1] and indices[lo] == s
    lo = ptr[t]
    hi = ptr[t + 1]
    while lo < hi:
        mid = (lo + hi) >> 1
        if idx[mid] < s:
            lo = mid + 1
        else:
            hi = mid
    if lo < ptr[t + 1] and idx[lo] == s:
        return adj, pid[lo]
    return adj, -1


@nb.njit(cache=True)
def _build_star(v, indptr, indices, eid, ptr, idx, pid, owner, perm_keys, T, Tc, inner):
    """Greedy star at centre v over pairs with owner == -1.  Returns (k, ni)."""
    d = indptr[v + 1] - indptr[v]
    o = np.argsort(perm_keys[indptr[v]:indptr[v + 1]])
    k = 0
    ni = 0
    for oj in range(d):
        p = indptr[v] + o[oj]
        if owner[eid[p]] != -1:
            continue
        t = indices[p]
        ok = True
        start = ni
        for j in range(k):
            adj, e2 = _find_pair(indptr, indices, ptr, idx, pid, t, T[j])
            if adj or e2 < 0 or owner[e2] != -1 or ni >= inner.shape[0]:
                ok = False
                break
            inner[ni] = e2
            ni += 1
        if not ok:
            ni = start
            continue
        T[k] = t
        Tc[k] = eid[p]
        k += 1
    return k, ni


@nb.njit(cache=True)
def _star_ls(n, indptr, indices, eid, ptr, idx, pid, npairs, init_order, iters, radius_cap,
             seed, pool_cap):
    np.random.seed(seed)
    owner = -np.ones(npairs, dtype=np.int64)
    maxd = 1
    for v in range(n):
        if indptr[v + 1] - indptr[v] > maxd:
            maxd = indptr[v + 1] - indptr[v]
    T = np.empty(maxd + 1, dtype=np.int64)
    Tc = np.empty(maxd + 1, dtype=np.int64)
    inner = np.empty(maxd * 8 + 64, dtype=np.int64)
    # star storage (append-only pool)
    smax = npairs // 3 + 16 + iters * 4
    s_center = np.empty(smax, dtype=np.int64)
    s_k = np.empty(smax, dtype=np.int64)
    s_start = np.empty(smax, dtype=np.int64)
    s_len = np.empty(smax, dtype=np.int64)
    s_alive = np.zeros(smax, dtype=np.bool_)
    pool = np.empty(pool_cap, dtype=np.int64)
    # stars per centre: linked lists
    head = -np.ones(n, dtype=np.int64)
    nxt = -np.ones(smax, dtype=np.int64)
    ns = 0
    npool = 0
    value = 0
    keys = np.random.random(indptr[n])

    def dummy():
        return 0

    for oi in range(n):
        v = init_order[oi]
        while True:
            k, ni = _build_star(v, indptr, indices, eid, ptr, idx, pid, owner, keys, T, Tc, inner)
            if k < 2 or ns >= smax or npool + k + ni > pool_cap:
                break
            s = ns
            ns += 1
            s_center[s] = v
            s_k[s] = k
            s_start[s] = npool
            s_len[s] = k + ni
            s_alive[s] = True
            for j in range(k):
                pool[npool] = Tc[j]
                owner[Tc[j]] = s
                npool += 1
            for j in range(ni):
                pool[npool] = inner[j]
                owner[inner[j]] = s
                npool += 1
            nxt[s] = head[v]
            head[v] = s
            value += k - 1
    best = value
    # local search
    region = np.empty(radius_cap, dtype=np.int64)
    inR = np.zeros(n, dtype=np.bool_)
    removed = np.empty(smax, dtype=np.int64)
    created = np.empty(smax, dtype=np.int64)
    for it in range(iters):
        if npool + 10 * maxd * maxd > pool_cap or ns + 4 * radius_cap > smax:
            # compact pool: rebuild from alive stars
            np2 = 0
            newpool = np.empty(pool_cap, dtype=np.int64)
            for s in range(ns):
                if s_alive[s]:
                    for q in range(s_len[s]):
                        newpool[np2 + q] = pool[s_start[s] + q]
                    s_start[s] = np2
                    np2 += s_len[s]
            pool[:np2] = newpool[:np2]
            npool = np2
            if ns + 4 * radius_cap > smax:
                break
        v0 = np.random.randint(n)
        if indptr[v0 + 1] - indptr[v0] == 0:
            continue
        # region: v0 and a random subset of its neighbours (+ their neighbours)
        nr = 1
        region[0] = v0
        inR[v0] = True
        h = 0
        while h < nr and nr < radius_cap:
            u = region[h]
            h += 1
            for p in range(indptr[u], indptr[u + 1]):
                w = indices[p]
                if not inR[w] and np.random.random() < 0.5:
                    inR[w] = True
                    region[nr] = w
                    nr += 1
                    if nr == radius_cap:
                        break
        # ruin: remove stars centred in the region
        nrem = 0
        old_val = 0
        for i in range(nr):
            u = region[i]
            s = head[u]
            while s != -1:
                if s_alive[s]:
                    s_alive[s] = False
                    removed[nrem] = s
                    nrem += 1
                    old_val += s_k[s] - 1
                    for q in range(s_len[s]):
                        owner[pool[s_start[s] + q]] = -1
                s = nxt[s]
            head[u] = -1
        # recreate in random order with fresh leaf keys
        for i in range(nr):
            u = region[i]
            for p in range(indptr[u], indptr[u + 1]):
                keys[p] = np.random.random()
        perm = np.random.permutation(nr)
        ncr = 0
        new_val = 0
        for i in range(nr):
            u = region[perm[i]]
            while True:
                k, ni = _build_star(u, indptr, indices, eid, ptr, idx, pid, owner, keys, T, Tc,
                                    inner)
                if k < 2 or ns >= smax or npool + k + ni > pool_cap:
                    break
                s = ns
                ns += 1
                s_center[s] = u
                s_k[s] = k
                s_start[s] = npool
                s_len[s] = k + ni
                s_alive[s] = True
                for j in range(k):
                    pool[npool] = Tc[j]
                    owner[Tc[j]] = s
                    npool += 1
                for j in range(ni):
                    pool[npool] = inner[j]
                    owner[inner[j]] = s
                    npool += 1
                nxt[s] = head[u]
                head[u] = s
                created[ncr] = s
                ncr += 1
                new_val += k - 1
        if new_val < old_val:
            # revert
            for c in range(ncr):
                s = created[c]
                s_alive[s] = False
                for q in range(s_len[s]):
                    owner[pool[s_start[s] + q]] = -1
            for i in range(nr):
                head[region[i]] = -1
            for c in range(nrem):
                s = removed[c]
                s_alive[s] = True
                for q in range(s_len[s]):
                    owner[pool[s_start[s] + q]] = s
                u = s_center[s]
                nxt[s] = head[u]
                head[u] = s
        else:
            value += new_val - old_val
        for i in range(nr):
            inR[region[i]] = False
    # export alive stars
    cnt = 0
    tot = 0
    for s in range(ns):
        if s_alive[s]:
            cnt += 1
            tot += s_len[s]
    rp = np.zeros(cnt + 1, dtype=np.int64)
    rpairs = np.empty(tot, dtype=np.int64)
    rk = np.empty(cnt, dtype=np.int64)
    c = 0
    z = 0
    for s in range(ns):
        if s_alive[s]:
            for q in range(s_len[s]):
                rpairs[z] = pool[s_start[s] + q]
                z += 1
            rk[c] = s_k[s]
            c += 1
            rp[c] = z
    return value, rp, rpairs, rk


def star_packing_ls(g: Graph, sup: Support | None = None, pgraph=None, iters: int = 200000,
                    region: int = 12, seed: int = 0):
    """Star packing improved by ruin-and-recreate local search (never decreases)."""
    from .lp import _pgraph
    sup = sup if sup is not None else build_support(g)
    if pgraph is None:
        pgraph = _pgraph(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id)
    ptr, idx, pid = pgraph
    rng = np.random.default_rng(seed)
    order = np.lexsort((rng.random(g.n), -g.degrees)).astype(np.int64)
    pool_cap = int(2 * sup.npairs + 10 * int(g.degrees.max()) ** 2 + 1000)
    return _star_ls(g.n, g.indptr, g.indices, sup.eid, ptr, idx, pid, sup.npairs, order,
                    iters, region, seed, pool_cap)


# P3-first packing: maximal packing of bad triangles, then star extensions that
# only use pairs no free bad triangle can use.

@nb.njit(cache=True)
def _try_p3(v, indptr, indices, eid, ptr, idx, pid, owner, keys, cap, T, Tc, inner):
    """Find a bad triangle centred at v on free pairs (first cap free neighbours in key
    order).  Returns (k, ni) = (2, 1) on success, (0, 0) otherwise."""
    d = indptr[v + 1] - indptr[v]
    o = np.argsort(keys[indptr[v]:indptr[v + 1]])
    fr = np.empty(min(d, cap), dtype=np.int64)
    nf = 0
    for oj in range(d):
        p = indptr[v] + o[oj]
        if owner[eid[p]] == -1:
            fr[nf] = p
            nf += 1
            if nf == fr.shape[0]:
                break
    for a in range(nf):
        pa = fr[a]
        ta = indices[pa]
        for b in range(a + 1, nf):
            pb = fr[b]
            tb = indices[pb]
            adj, e2 = _find_pair(indptr, indices, ptr, idx, pid, ta, tb)
            if adj or e2 < 0 or owner[e2] != -1:
                continue
            T[0] = ta
            T[1] = tb
            Tc[0] = eid[pa]
            Tc[1] = eid[pb]
            inner[0] = e2
            return 2, 1
    return 0, 0


@nb.njit(cache=True)
def _try_extend(v, k, ni, indptr, indices, eid, ptr, idx, pid, owner, keys, T, Tc, inner):
    """Greedily add leaves to the star (v; T[:k]) using free pairs only."""
    d = indptr[v + 1] - indptr[v]
    o = np.argsort(keys[indptr[v]:indptr[v + 1]])
    for oj in range(d):
        p = indptr[v] + o[oj]
        if owner[eid[p]] != -1:
            continue
        t = indices[p]
        ok = True
        start = ni
        for j in range(k):
            if T[j] == t:
                ok = False
                break
            adj, e2 = _find_pair(indptr, indices, ptr, idx, pid, t, T[j])
            if adj or e2 < 0 or owner[e2] != -1 or ni >= inner.shape[0]:
                ok = False
                break
            inner[ni] = e2
            ni += 1
        if not ok:
            ni = start
            continue
        if k >= T.shape[0]:
            break
        T[k] = t
        Tc[k] = eid[p]
        k += 1
    return k, ni


@nb.njit(cache=True)
def _p3x_ls(n, indptr, indices, eid, ptr, idx, pid, npairs, iters, radius_cap, seed, cap,
            pool_cap):
    """P3-first star packing with ruin-and-recreate local search.

    Construction (globally and inside every recreated region): (1) a maximal
    packing of bad triangles on free pairs, (2) star extensions that add a
    leaf using only free pairs.  After (1) no bad triangle has three free
    pairs, so an extension never blocks a triangle.  A local-search step
    removes all stars centred in a random region and rebuilds it; it is kept
    iff the value does not decrease."""
    np.random.seed(seed)
    owner = -np.ones(npairs, dtype=np.int64)
    maxd = 2
    for v in range(n):
        if indptr[v + 1] - indptr[v] > maxd:
            maxd = indptr[v + 1] - indptr[v]
    keys = np.random.random(indptr[n])
    T = np.empty(maxd + 1, dtype=np.int64)
    Tc = np.empty(maxd + 1, dtype=np.int64)
    inner = np.empty(min(maxd * (maxd + 1) // 2 + 8, 5000000), dtype=np.int64)
    smax = npairs // 3 + 8 + 4 * (radius_cap * (maxd + 1) + 16)
    s_center = np.empty(smax, dtype=np.int64)
    s_k = np.zeros(smax, dtype=np.int64)
    s_ps = np.zeros(smax, dtype=np.int64)   # pair segment start
    s_pl = np.zeros(smax, dtype=np.int64)   # pair segment length
    s_ls = np.zeros(smax, dtype=np.int64)   # leaf segment start
    s_alive = np.zeros(smax, dtype=np.bool_)
    ppool = np.empty(pool_cap, dtype=np.int64)
    lpool = np.empty(pool_cap, dtype=np.int64)
    npp = 0
    nlp = 0
    ns = 0
    head = -np.ones(n, dtype=np.int64)
    nxt = -np.ones(smax, dtype=np.int64)
    value = 0
    region = np.empty(max(radius_cap, n), dtype=np.int64)
    inR = np.zeros(n, dtype=np.bool_)
    removed = np.empty(smax, dtype=np.int64)
    created = np.empty(smax, dtype=np.int64)
    old_k = np.empty(smax, dtype=np.int64)

    for phase_it in range(iters + 1):
        # ---------------- choose region ----------------
        if phase_it == 0:
            nr = n
            for i in range(n):
                region[i] = i
        else:
            v0 = np.random.randint(n)
            if indptr[v0 + 1] - indptr[v0] == 0:
                continue
            nr = 1
            region[0] = v0
            inR[v0] = True
            h = 0
            while h < nr and nr < radius_cap:
                u = region[h]
                h += 1
                for p in range(indptr[u], indptr[u + 1]):
                    w = indices[p]
                    if not inR[w] and np.random.random() < 0.5:
                        inR[w] = True
                        region[nr] = w
                        nr += 1
                        if nr == radius_cap:
                            break
            for i in range(nr):
                inR[region[i]] = False
        # compaction of star ids, lists and pools if needed
        if phase_it > 0 and (ns + nr * maxd + 8 > smax or
                             npp + nr * (maxd + 8) * 4 + 64 > pool_cap or
                             nlp + nr * (maxd + 4) * 2 + 64 > pool_cap):
            newid = -np.ones(ns, dtype=np.int64)
            c2 = 0
            for s in range(ns):
                if s_alive[s]:
                    newid[s] = c2
                    s_center[c2] = s_center[s]
                    s_k[c2] = s_k[s]
                    s_ps[c2] = s_ps[s]
                    s_pl[c2] = s_pl[s]
                    s_ls[c2] = s_ls[s]
                    s_alive[c2] = True
                    c2 += 1
            for s in range(c2, ns):
                s_alive[s] = False
            for e in range(npairs):
                if owner[e] >= 0:
                    owner[e] = newid[owner[e]]
            ns = c2
            for v in range(n):
                head[v] = -1
            for s in range(ns):
                u = s_center[s]
                nxt[s] = head[u]
                head[u] = s
            if ns + nr * maxd + 8 > smax:
                break
            np2 = 0
            nl2 = 0
            tmpp = ppool.copy()
            tmpl = lpool.copy()
            for s in range(ns):
                if s_alive[s]:
                    for q in range(s_pl[s]):
                        ppool[np2 + q] = tmpp[s_ps[s] + q]
                    s_ps[s] = np2
                    np2 += s_pl[s]
                    for q in range(s_k[s]):
                        lpool[nl2 + q] = tmpl[s_ls[s] + q]
                    s_ls[s] = nl2
                    nl2 += s_k[s]
            npp = np2
            nlp = nl2
            if npp + nr * (maxd + 8) * 4 + 64 > pool_cap:
                break
        # ---------------- ruin ----------------
        nrem = 0
        old_val = 0
        if phase_it > 0:
            for i in range(nr):
                u = region[i]
                s = head[u]
                while s != -1:
                    if s_alive[s]:
                        s_alive[s] = False
                        removed[nrem] = s
                        old_k[nrem] = s_k[s]
                        nrem += 1
                        old_val += s_k[s] - 1
                        for q in range(s_pl[s]):
                            owner[ppool[s_ps[s] + q]] = -1
                    s = nxt[s]
                head[u] = -1
            for i in range(nr):
                u = region[i]
                for p in range(indptr[u], indptr[u + 1]):
                    keys[p] = np.random.random()
        # ---------------- recreate: triangles ----------------
        ncr = 0
        perm = np.random.permutation(nr)
        for i in range(nr):
            u = region[perm[i]]
            while True:
                k, ni = _try_p3(u, indptr, indices, eid, ptr, idx, pid, owner, keys, cap, T, Tc,
                                inner)
                if k == 0:
                    break
                if ns >= smax:
                    break
                s = ns
                ns += 1
                s_center[s] = u
                s_k[s] = 2
                s_ps[s] = npp
                s_pl[s] = 3
                ppool[npp] = Tc[0]
                ppool[npp + 1] = Tc[1]
                ppool[npp + 2] = inner[0]
                npp += 3
                s_ls[s] = nlp
                lpool[nlp] = T[0]
                lpool[nlp + 1] = T[1]
                nlp += 2
                s_alive[s] = True
                owner[Tc[0]] = s
                owner[Tc[1]] = s
                owner[inner[0]] = s
                nxt[s] = head[u]
                head[u] = s
                created[ncr] = s
                ncr += 1
        # ---------------- recreate: extensions ----------------
        new_val = 0
        for c in range(ncr):
            s = created[c]
            u = s_center[s]
            k = s_k[s]
            for j in range(k):
                T[j] = lpool[s_ls[s] + j]
            k2, ni2 = _try_extend(u, k, 0, indptr, indices, eid, ptr, idx, pid, owner, keys,
                                  T, Tc, inner)
            if k2 > k:
                # re-store the star with its new leaves and pairs
                newps = npp
                for q in range(s_pl[s]):
                    ppool[npp] = ppool[s_ps[s] + q]
                    npp += 1
                for j in range(k, k2):
                    ppool[npp] = Tc[j]
                    owner[Tc[j]] = s
                    npp += 1
                for q in range(ni2):
                    ppool[npp] = inner[q]
                    owner[inner[q]] = s
                    npp += 1
                s_pl[s] = npp - newps
                s_ps[s] = newps
                newls = nlp
                for j in range(k2):
                    lpool[nlp] = T[j]
                    nlp += 1
                s_ls[s] = newls
                s_k[s] = k2
            new_val += s_k[s] - 1
        # ---------------- accept / revert ----------------
        if phase_it == 0 or new_val >= old_val:
            value += new_val - old_val
        else:
            for c in range(ncr):
                s = created[c]
                s_alive[s] = False
                for q in range(s_pl[s]):
                    owner[ppool[s_ps[s] + q]] = -1
            for i in range(nr):
                head[region[i]] = -1
            for c in range(nrem):
                s = removed[c]
                s_alive[s] = True
                for q in range(s_pl[s]):
                    owner[ppool[s_ps[s] + q]] = s
                u = s_center[s]
                nxt[s] = head[u]
                head[u] = s
        # rebuild head lists lazily: dead stars are skipped when walking
    # export
    cnt = 0
    tot = 0
    for s in range(ns):
        if s_alive[s]:
            cnt += 1
            tot += s_pl[s]
    rp = np.zeros(cnt + 1, dtype=np.int64)
    rpairs = np.empty(tot, dtype=np.int64)
    rk = np.empty(cnt, dtype=np.int64)
    c = 0
    z = 0
    for s in range(ns):
        if s_alive[s]:
            # centre pairs first: the first two pool entries are centre pairs, then the
            # triangle's leaf pair, then per extension: centre pair + leaf pairs.
            # Reorder: centre pairs (+1) then leaf pairs (-1).
            k = s_k[s]
            base = s_ps[s]
            # centre pairs are those incident to the centre
            for q in range(s_pl[s]):
                e = ppool[base + q]
                rpairs[z] = e
                z += 1
            rk[c] = k
            c += 1
            rp[c] = z
    return value, rp, rpairs, rk


def star_packing_p3x(g: Graph, sup: Support | None = None, pgraph=None, iters: int = 200000,
                     region: int = 12, seed: int = 0, cap: int = 64):
    """P3-first star packing with ruin-and-recreate local search (never decreases).

    Returns (value, row_ptr, row_pairs, row_k) like star_packing_ls."""
    from .lp import _pgraph
    sup = sup if sup is not None else build_support(g)
    if pgraph is None:
        pgraph = _pgraph(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id)
    ptr, idx, pid = pgraph
    pool_cap = int(3 * sup.npairs + 40 * (int(g.degrees.max()) + 8) * region + 1000)
    return _p3x_ls(g.n, g.indptr, g.indices, sup.eid, ptr, idx, pid, sup.npairs, iters,
                   region, seed, cap, pool_cap)
