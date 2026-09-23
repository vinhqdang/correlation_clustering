"""Sparse LP relaxation, LP rounding and exact ILP on the distance-2 support.

Variables x_e in [0, 1] (distance / "separated" indicator) for e in P = E+ u N2;
far pairs (G+-distance >= 3) are fixed to x = 1.  By the diameter-2 lemma this
keeps every optimal clustering feasible, so the LP value is a lower bound on
OPT, and because all triangle inequalities involving far pairs are enforced
the solution (extended by x = 1 on far pairs) is a feasible solution of the
standard metric LP.  Hence every LP-rounding guarantee proved against the
standard LP applies verbatim.

Triangle inequalities are generated lazily (cutting planes).
"""
from __future__ import annotations

import time

import highspy
import numpy as np
import numba as nb

from .graph import Graph
from .support import Support, build_support, _find, wedge_triangles


@nb.njit(cache=True)
def _pgraph(n, indptr, indices, eid, n2ptr, n2idx, n2id):
    """CSR of the support graph P (union of E+ and N2 rows, sorted)."""
    ptr = np.zeros(n + 1, dtype=np.int64)
    for u in range(n):
        ptr[u + 1] = ptr[u] + (indptr[u + 1] - indptr[u]) + (n2ptr[u + 1] - n2ptr[u])
    idx = np.empty(ptr[n], dtype=np.int32)
    pid = np.empty(ptr[n], dtype=np.int64)
    for u in range(n):
        i = indptr[u]
        j = n2ptr[u]
        k = ptr[u]
        while i < indptr[u + 1] or j < n2ptr[u + 1]:
            if j >= n2ptr[u + 1] or (i < indptr[u + 1] and indices[i] < n2idx[j]):
                idx[k] = indices[i]
                pid[k] = eid[i]
                i += 1
            else:
                idx[k] = n2idx[j]
                pid[k] = n2id[j]
                j += 1
            k += 1
    return ptr, idx, pid


@nb.njit(cache=True)
def _separate(n, ptr, idx, pid, x, eps, cap):
    """Find violated triangle inequalities.

    Returns rows as (type, e1, e2, e3, viol):
      type 0:  x[e3] - x[e1] - x[e2] <= 0   (e3 = uv, e1 = uw, e2 = wv)
      type 1:  x[e1] + x[e2] >= 1           (uv far)
    """
    out_t = np.empty(cap, dtype=np.int8)
    out_a = np.empty(cap, dtype=np.int64)
    out_b = np.empty(cap, dtype=np.int64)
    out_c = np.empty(cap, dtype=np.int64)
    out_v = np.empty(cap, dtype=np.float64)
    k = 0
    close = np.empty(n, dtype=np.int64)
    closep = np.empty(n, dtype=np.int64)
    for w in range(n):
        nc = 0
        for p in range(ptr[w], ptr[w + 1]):
            if x[pid[p]] < 1.0 - eps:
                close[nc] = idx[p]
                closep[nc] = pid[p]
                nc += 1
        for i in range(nc):
            u = close[i]
            xu = x[closep[i]]
            for j in range(i + 1, nc):
                v = close[j]
                s = xu + x[closep[j]]
                if s >= 1.0 - eps:
                    continue
                r = _find(ptr, idx, u, v)
                if r >= 0:
                    viol = x[pid[r]] - s
                    typ = 0
                    e3 = pid[r]
                else:
                    viol = 1.0 - s
                    typ = 1
                    e3 = -1
                if viol > eps:
                    if k < cap:
                        out_t[k] = typ
                        out_a[k] = closep[i]
                        out_b[k] = closep[j]
                        out_c[k] = e3
                        out_v[k] = viol
                    k += 1
    kk = min(k, cap)
    return out_t[:kk], out_a[:kk], out_b[:kk], out_c[:kk], out_v[:kk], k


class SparseLP:
    """Cutting-plane solver for the metric LP on the distance-2 support."""

    def __init__(self, g: Graph, sup: Support | None = None, integral: bool = False,
                 verbose: bool = False, threads: int = 1):
        self.g = g
        self.sup = sup if sup is not None else build_support(g)
        s = self.sup
        self.ptr, self.idx, self.pid = _pgraph(g.n, g.indptr, g.indices, s.eid,
                                               s.n2ptr, s.n2idx, s.n2id)
        self.integral = integral
        h = highspy.Highs()
        h.setOptionValue("output_flag", bool(verbose))
        h.setOptionValue("threads", int(threads))
        N = s.npairs
        cost = np.concatenate([np.ones(s.npos), -np.ones(s.nneg)])
        h.addVars(N, np.zeros(N), np.ones(N))
        h.changeColsCost(N, np.arange(N, dtype=np.int32), cost)
        if integral:
            h.changeColsIntegrality(N, np.arange(N, dtype=np.int32),
                                    np.array([highspy.HighsVarType.kInteger] * N))
        self.h = h
        self.const = float(s.nneg)
        self.nrows = 0
        self.rounds = 0
        # initial rows: every bad-triangle inequality x_uv <= x_uw + x_wv
        ta, tb, tc = wedge_triangles(s)
        bad = tc >= s.npos
        self._add(np.zeros(int(bad.sum()), dtype=np.int8), ta[bad], tb[bad], tc[bad])

    def _add(self, typ, a, b, c):
        k = len(typ)
        if k == 0:
            return
        is0 = typ == 0
        nnz_per = np.where(is0, 3, 2)
        starts = np.zeros(k, dtype=np.int32)
        starts[1:] = np.cumsum(nnz_per)[:-1]
        nnz = int(nnz_per.sum())
        index = np.empty(nnz, dtype=np.int32)
        value = np.empty(nnz)
        pos = starts
        index[pos] = a
        value[pos] = 1.0
        index[pos + 1] = b
        value[pos + 1] = 1.0
        i0 = np.flatnonzero(is0)
        index[pos[i0] + 2] = c[i0]
        value[pos[i0] + 2] = -1.0
        lower = np.where(is0, 0.0, 1.0)
        upper = np.full(k, highspy.kHighsInf)
        self.h.addRows(k, lower, upper, nnz, starts, index, value)
        self.nrows += k

    def solve(self, eps: float = 1e-6, max_rounds: int = 200, cap: int = 200000,
              time_limit: float = 3600.0):
        t0 = time.time()
        while True:
            self.h.setOptionValue("time_limit", max(1.0, time_limit - (time.time() - t0)))
            self.h.run()
            x = np.asarray(self.h.getSolution().col_value)
            self.rounds += 1
            typ, a, b, c, v, found = _separate(self.g.n, self.ptr, self.idx, self.pid,
                                               x, eps, cap)
            if found == 0 or self.rounds >= max_rounds or time.time() - t0 > time_limit:
                break
            if found > cap:
                order = np.argsort(-v)
                typ, a, b, c = typ[order], a[order], b[order], c[order]
            self._add(typ, a, b, c)
        self.x = np.clip(x, 0.0, 1.0)
        self.value = self.const + float(self.h.getInfo().objective_function_value)
        self.converged = found == 0
        self.time = time.time() - t0
        return self.value

    def dual_bound(self) -> float:
        """For the MIP: best proven lower bound; for the LP: the LP value."""
        if self.integral:
            return self.const + float(self.h.getInfo().mip_dual_bound)
        return self.value

    def clustering(self) -> np.ndarray:
        """Clusters from an integral solution (components of x_e = 0 pairs)."""
        from .graph import from_edges, components
        s = self.sup
        join = self.x < 0.5
        e = np.stack([s.pu[join], s.pv[join]], axis=1)
        return components(from_edges(self.g.n, e))[0]


# --------------------------------------------------------------------------
# LP rounding
# --------------------------------------------------------------------------

CMSY_A = 0.19
CMSY_B = 0.5095


@nb.njit(cache=True)
def _round(n, ptr, idx, pid, x, npos, order, rnd, scheme, a, b):
    labels = np.full(n, -1, dtype=np.int64)
    c = 0
    k = 0
    for i in range(n):
        u = order[i]
        if labels[u] >= 0:
            continue
        labels[u] = c
        for p in range(ptr[u], ptr[u + 1]):
            v = idx[p]
            if labels[v] >= 0:
                continue
            e = pid[p]
            xe = x[e]
            if scheme == 0 or e >= npos:
                f = xe
            else:
                if xe < a:
                    f = 0.0
                elif xe < b:
                    f = ((xe - a) / (b - a)) ** 2
                else:
                    f = 1.0
            if rnd[k % rnd.shape[0]] < 1.0 - f:
                labels[v] = c
            k += 1
        c += 1
    return labels


def lp_pivot(g: Graph, sup: Support, x: np.ndarray, scheme: str = "cmsy", rng=None,
             pgraph=None) -> np.ndarray:
    """Pivot-based rounding of a metric x on P (far pairs have x = 1).

    scheme="acn":  join with probability 1 - x   (Ailon-Charikar-Newman, 2.5)
    scheme="cmsy": positive pairs use f+(x)=((x-a)/(b-a))^2 on [a,b], negative
                   pairs f-(x)=x (Chawla-Makarychev-Schramm-Yaroslavtsev)."""
    rng = np.random.default_rng(rng)
    if pgraph is None:
        pgraph = _pgraph(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id)
    ptr, idx, pid = pgraph
    order = rng.permutation(g.n).astype(np.int64)
    rnd = rng.random(max(1, int(ptr[-1])))
    return _round(g.n, ptr, idx, pid, np.asarray(x, dtype=np.float64), sup.npos, order,
                  rnd, 0 if scheme == "acn" else 1, CMSY_A, CMSY_B)


def lp_cost(sup: Support, x: np.ndarray) -> float:
    """Fractional cost of x (far pairs contribute 0)."""
    return float(x[: sup.npos].sum() + (1.0 - x[sup.npos:]).sum())


# --------------------------------------------------------------------------
# metric repair
# --------------------------------------------------------------------------

@nb.njit(cache=True)
def _repair(n, ptr, idx, pid, x, npairs):
    """d_uv = min(1, shortest-path distance w.r.t. x) on the pairs of P.

    d is a (pseudo)metric bounded by 1.  Pairs outside P keep d = 1, which is
    only consistent if no short path exists; we therefore also return the
    total "far deficit" sum over far pairs of (1 - d_uv)_+ restricted to far
    pairs reachable within distance < 1 (so the caller can account for it)."""
    d = x.copy()
    dist = np.full(n, 2.0)
    heap_v = np.empty(ptr[n] + n + 1, dtype=np.int64)
    heap_d = np.empty(ptr[n] + n + 1, dtype=np.float64)
    seen = np.empty(n, dtype=np.int64)
    far_deficit = 0.0
    far_pairs = 0
    for s in range(n):
        ns = 0
        hs = 0
        heap_v[0] = s
        heap_d[0] = 0.0
        hs = 1
        dist[s] = 0.0
        seen[ns] = s
        ns += 1
        while hs > 0:
            # pop min
            u = heap_v[0]
            du = heap_d[0]
            hs -= 1
            if hs > 0:
                lv = heap_v[hs]
                ld = heap_d[hs]
                i = 0
                while True:
                    l = 2 * i + 1
                    if l >= hs:
                        break
                    r = l + 1
                    c = l
                    if r < hs and heap_d[r] < heap_d[l]:
                        c = r
                    if heap_d[c] < ld:
                        heap_v[i] = heap_v[c]
                        heap_d[i] = heap_d[c]
                        i = c
                    else:
                        break
                heap_v[i] = lv
                heap_d[i] = ld
            if du > dist[u]:
                continue
            for p in range(ptr[u], ptr[u + 1]):
                v = idx[p]
                nd = du + x[pid[p]]
                if nd < 1.0 - 1e-9 and nd < dist[v] - 1e-12:
                    if dist[v] >= 2.0:
                        seen[ns] = v
                        ns += 1
                    dist[v] = nd
                    # push
                    i = hs
                    hs += 1
                    while i > 0:
                        par = (i - 1) >> 1
                        if heap_d[par] > nd:
                            heap_v[i] = heap_v[par]
                            heap_d[i] = heap_d[par]
                            i = par
                        else:
                            break
                    heap_v[i] = v
                    heap_d[i] = nd
        # update pairs of P incident to s, count far pairs
        for p in range(ptr[s], ptr[s + 1]):
            v = idx[p]
            if dist[v] < d[pid[p]]:
                d[pid[p]] = dist[v]
        for i in range(ns):
            v = seen[i]
            if v > s:
                if _find(ptr, idx, s, v) < 0:
                    far_deficit += 1.0 - dist[v]
                    far_pairs += 1
        for i in range(ns):
            dist[seen[i]] = 2.0
    return d, far_deficit, far_pairs


def repair_metric(g: Graph, sup: Support, x: np.ndarray, pgraph=None):
    """Project an arbitrary x on P to a feasible LP point d <= x.

    Returns (d, far_deficit, n_far_close): d restricted to P, and the extra
    fractional cost 1 - d_uv of far pairs that became closer than 1."""
    if pgraph is None:
        pgraph = _pgraph(g.n, g.indptr, g.indices, sup.eid, sup.n2ptr, sup.n2idx, sup.n2id)
    ptr, idx, pid = pgraph
    return _repair(g.n, ptr, idx, pid, np.clip(np.asarray(x, dtype=np.float64), 0, 1),
                   sup.npairs)


# --------------------------------------------------------------------------
# star (2-partition) inequalities
# --------------------------------------------------------------------------
#
# For a vertex v and a set T of other vertices, every clustering satisfies
#     sum_{t in T} y_vt - sum_{t < t' in T} y_tt' <= 1,     y = 1 - x
# (if v's cluster contains k vertices of T the left side is k - k(k-1)/2 <= 1).
# With |T| = 2 these are triangle inequalities.  Pairs outside P have y = 0 in
# every optimal solution and drop out.  In x-variables:
#   -sum_t x_vt + sum_{tt' in P} x_tt' <= 1 - |T| + |{tt' in P}|.

@nb.njit(cache=True)
def _separate_stars(n, ptr, idx, pid, x, eps, max_t, cap_rows, cap_nnz):
    rows_ptr = np.zeros(cap_rows + 1, dtype=np.int64)
    rows_idx = np.empty(cap_nnz, dtype=np.int64)
    rows_val = np.empty(cap_nnz, dtype=np.float64)
    rows_ub = np.empty(cap_rows, dtype=np.float64)
    rows_viol = np.empty(cap_rows, dtype=np.float64)
    nr = 0
    nz = 0
    cand = np.empty(n, dtype=np.int64)
    candy = np.empty(n, dtype=np.float64)
    candp = np.empty(n, dtype=np.int64)
    T = np.empty(max_t, dtype=np.int64)
    Tp = np.empty(max_t, dtype=np.int64)
    for v in range(n):
        k = 0
        for p in range(ptr[v], ptr[v + 1]):
            y = 1.0 - x[pid[p]]
            if y > eps:
                cand[k] = idx[p]
                candy[k] = y
                candp[k] = pid[p]
                k += 1
        if k < 3:
            continue
        order = np.argsort(-candy[:k])
        nt = 0
        lhs = 0.0
        for oi in range(k):
            i = order[oi]
            t = cand[i]
            gain = candy[i]
            for j in range(nt):
                r = _find(ptr, idx, t, T[j])
                if r >= 0:
                    gain -= 1.0 - x[pid[r]]
            if gain > eps:
                T[nt] = t
                Tp[nt] = candp[i]
                nt += 1
                lhs += gain
                if nt == max_t:
                    break
        if nt >= 3 and lhs > 1.0 + eps:
            # build row
            need = nt + nt * (nt - 1) // 2
            if nr >= cap_rows or nz + need > cap_nnz:
                break
            ub = 1.0 - nt
            for j in range(nt):
                rows_idx[nz] = Tp[j]
                rows_val[nz] = -1.0
                nz += 1
            for a in range(nt):
                for b in range(a + 1, nt):
                    r = _find(ptr, idx, T[a], T[b])
                    if r >= 0:
                        rows_idx[nz] = pid[r]
                        rows_val[nz] = 1.0
                        nz += 1
                        ub += 1.0
            rows_ub[nr] = ub
            rows_viol[nr] = lhs - 1.0
            nr += 1
            rows_ptr[nr] = nz
    return rows_ptr[:nr + 1], rows_idx[:nz], rows_val[:nz], rows_ub[:nr], rows_viol[:nr]


def _add_general(self, rptr, ridx, rval, lower, upper):
    k = len(lower)
    if k == 0:
        return
    self.h.addRows(k, np.asarray(lower, dtype=np.float64), np.asarray(upper, dtype=np.float64),
                   int(rptr[-1]), rptr[:-1].astype(np.int32), ridx.astype(np.int32),
                   rval.astype(np.float64))
    self.nrows += k


def _solve_with_stars(self, eps: float = 1e-6, max_rounds: int = 200, time_limit: float = 3600.0,
                      max_t: int = 50, star_cap: int = 50000, cap: int = 200000,
                      verbose: bool = False):
    """Cutting planes with triangle and star inequalities."""
    t0 = time.time()
    self.history = []
    while True:
        self.h.setOptionValue("time_limit", max(1.0, time_limit - (time.time() - t0)))
        self.h.run()
        x = np.asarray(self.h.getSolution().col_value)
        val = self.const + float(self.h.getInfo().objective_function_value)
        self.rounds += 1
        self.history.append((self.rounds, time.time() - t0, val, self.nrows))
        if verbose:
            print(self.rounds, round(time.time() - t0, 1), val, self.nrows, flush=True)
        typ, a, b, c, v, found = _separate(self.g.n, self.ptr, self.idx, self.pid, x, eps, cap)
        if found > cap:
            o = np.argsort(-v)
            typ, a, b, c = typ[o], a[o], b[o], c[o]
        self._add(typ, a, b, c)
        nstar = 0
        if found == 0:
            rp, ri, rv, ub, viol = _separate_stars(self.g.n, self.ptr, self.idx, self.pid, x,
                                                   1e-4, max_t, star_cap, star_cap * 60)
            nstar = len(ub)
            _add_general(self, rp, ri, rv, np.full(nstar, -highspy.kHighsInf), ub)
        if (found == 0 and nstar == 0) or self.rounds >= max_rounds or \
                time.time() - t0 > time_limit:
            break
    self.x = np.clip(x, 0.0, 1.0)
    self.value = val
    self.converged = found == 0 and nstar == 0
    self.time = time.time() - t0
    return self.value


SparseLP.solve_with_stars = _solve_with_stars
SparseLP.add_general = _add_general
