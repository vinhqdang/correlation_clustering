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
