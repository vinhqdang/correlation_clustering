"""Scalable certified lower bounds by dual block-coordinate ascent.

Write the relaxation over the distance-2 support P as

    min  K + sum_e c_e x_e   s.t.  A x >= b,  0 <= x <= 1,

with c_e = +1 on E+, -1 on N2 and K = |N2|; the rows are triangle
inequalities (including the ones through far pairs) and star inequalities.
For ANY multipliers y >= 0 on ANY set of valid rows, weak duality gives

    LB(y) = K + b^T y + sum_e min(0, r_e),     r = c - A^T y,         (*)

because x ranges over the box.  (*) only needs the running quantities
``r`` (one number per pair) and ``b^T y``; rows never have to be stored.

A block step picks a vertex set B, takes the current reduced costs r on the
pairs of P inside B as costs, and solves the LP restricted to those pairs
with lazily separated triangle and star rows.  Its optimal duals y_B are
added to y, which increases (*) by

    LP_B(r) - sum_{e in P_B} min(0, r_e)  >= 0 .

The bound is therefore monotone and valid after every block, blocks of one
partition are independent, and successive sweeps over different partitions
let rows cross earlier block boundaries.
"""
from __future__ import annotations

import time

import highspy
import numpy as np
import numba as nb

from .graph import Graph
from .support import Support, build_support
from .lp import _pgraph, _separate, _separate_stars, _separate_subgraphs


@nb.njit(cache=True)
def _local_pgraph(nodes, pos, ptr, idx, pid):
    """Pair graph of P restricted to the vertex set ``nodes`` (pos = local index or -1)."""
    k = nodes.shape[0]
    lptr = np.zeros(k + 1, dtype=np.int64)
    for i in range(k):
        v = nodes[i]
        c = 0
        for p in range(ptr[v], ptr[v + 1]):
            if pos[idx[p]] >= 0:
                c += 1
        lptr[i + 1] = lptr[i] + c
    lidx = np.empty(lptr[k], dtype=np.int32)
    lgid = np.empty(lptr[k], dtype=np.int64)
    for i in range(k):
        v = nodes[i]
        q = lptr[i]
        for p in range(ptr[v], ptr[v + 1]):
            j = pos[idx[p]]
            if j >= 0:
                lidx[q] = j
                lgid[q] = pid[p]
                q += 1
        # rows are sorted by global id order of neighbours; local ids must be sorted
        seg = lidx[lptr[i]:lptr[i + 1]]
        o = np.argsort(seg)
        lidx[lptr[i]:lptr[i + 1]] = seg[o]
        g2 = lgid[lptr[i]:lptr[i + 1]].copy()
        lgid[lptr[i]:lptr[i + 1]] = g2[o]
    # local pair ids: compress global ids appearing in the block
    return lptr, lidx, lgid


class _BlockLP:
    def __init__(self, costs, n, ptr, idx, pid, npos=None):
        self.n, self.ptr, self.idx, self.pid = n, ptr, idx, pid
        self.npos = npos
        N = costs.shape[0]
        h = highspy.Highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("threads", 1)
        h.addVars(N, np.zeros(N), np.ones(N))
        h.changeColsCost(N, np.arange(N, dtype=np.int32), costs)
        self.h = h
        self.costs = costs
        # rows kept in >= form for dual bookkeeping
        self.r_ptr = [0]
        self.r_idx = []
        self.r_val = []
        self.r_b = []

    def _push(self, starts, index, value, lower):
        k = len(lower)
        if k == 0:
            return
        self.h.addRows(k, lower, np.full(k, highspy.kHighsInf), len(index),
                       starts.astype(np.int32), index.astype(np.int32), value)
        off = self.r_ptr[-1]
        ends = np.append(starts[1:], len(index))
        self.r_ptr.extend((off + ends).tolist())
        self.r_idx.append(index.astype(np.int64))
        self.r_val.append(value.astype(np.float64))
        self.r_b.append(np.asarray(lower, dtype=np.float64))

    def add_triangles(self, typ, a, b, c):
        k = len(typ)
        if k == 0:
            return
        is0 = typ == 0
        per = np.where(is0, 3, 2)
        starts = np.zeros(k, dtype=np.int64)
        starts[1:] = np.cumsum(per)[:-1]
        nnz = int(per.sum())
        index = np.empty(nnz, dtype=np.int64)
        value = np.empty(nnz)
        index[starts] = a
        value[starts] = 1.0
        index[starts + 1] = b
        value[starts + 1] = 1.0
        i0 = np.flatnonzero(is0)
        index[starts[i0] + 2] = c[i0]
        value[starts[i0] + 2] = -1.0
        self._push(starts, index, value, np.where(is0, 0.0, 1.0))

    def add_stars(self, rp, ri, rv, ub):
        # -sum x_vt + sum x_tt' <= ub   ->   sum x_vt - sum x_tt' >= -ub
        if len(ub) == 0:
            return
        self._push(rp[:-1].astype(np.int64), ri.astype(np.int64), -rv, -np.asarray(ub))

    def add_subgraphs(self, rp, ri, rv, lb):
        if len(lb) == 0:
            return
        self._push(rp[:-1].astype(np.int64), ri.astype(np.int64), rv, np.asarray(lb))

    def solve(self, eps=1e-6, max_rounds=100, time_limit=600.0, stars=True, max_t=50,
              star_cap=20000, cap=100000, subgraphs=False, hmax=8, sub_cap=20000, seed=0):
        t0 = time.time()
        rounds = 0
        while True:
            self.h.setOptionValue("time_limit", max(1.0, time_limit - (time.time() - t0)))
            self.h.run()
            x = np.asarray(self.h.getSolution().col_value)
            rounds += 1
            typ, a, b, c, v, found = _separate(self.n, self.ptr, self.idx, self.pid, x, eps, cap)
            if found > cap:
                o = np.argsort(-v)
                typ, a, b, c = typ[o], a[o], b[o], c[o]
            self.add_triangles(typ, a, b, c)
            ns = 0
            if found == 0 and stars:
                rp, ri, rv, ub, _ = _separate_stars(self.n, self.ptr, self.idx, self.pid, x,
                                                    1e-4, max_t, star_cap, star_cap * 60)
                ns = len(ub)
                self.add_stars(rp, ri, rv, ub)
            if found == 0 and ns == 0 and subgraphs and self.npos is not None:
                order = np.random.default_rng(seed + rounds).permutation(self.n).astype(np.int64)
                rp, ri, rv, lbs = _separate_subgraphs(self.n, self.ptr, self.idx, self.pid, x,
                                                      self.npos, 1e-4, hmax, sub_cap, order, True)
                ns = len(lbs)
                self.add_subgraphs(rp, ri, rv, lbs)
            if (found == 0 and ns == 0) or rounds >= max_rounds or \
                    time.time() - t0 > time_limit:
                break
        self.rounds = rounds
        self.converged = found == 0 and ns == 0
        return x

    def duals(self):
        y = np.asarray(self.h.getSolution().row_dual, dtype=np.float64)
        y = np.maximum(y, 0.0)
        return y

    def rows(self):
        if not self.r_idx:
            return (np.zeros(1, dtype=np.int64), np.zeros(0, dtype=np.int64),
                    np.zeros(0), np.zeros(0))
        return (np.asarray(self.r_ptr, dtype=np.int64), np.concatenate(self.r_idx),
                np.concatenate(self.r_val), np.concatenate(self.r_b))


def _gather(ptr, sel):
    """Flat positions of the rows ``sel`` of a CSR and the row index of each."""
    lens = ptr[sel + 1] - ptr[sel]
    rowof = np.repeat(np.arange(len(sel)), lens)
    starts = np.repeat(ptr[sel], lens)
    within = np.arange(int(lens.sum())) - np.repeat(np.cumsum(lens) - lens, lens)
    return starts + within, rowof


@nb.njit(cache=True)
def _apply_duals(rptr, ridx, rval, y, r):
    """r -= A^T y (in place, local ids); returns b^T y is computed outside."""
    for i in range(rptr.shape[0] - 1):
        yi = y[i]
        if yi == 0.0:
            continue
        for p in range(rptr[i], rptr[i + 1]):
            r[ridx[p]] -= rval[p] * yi


class BlockDualBound:
    """Anytime-valid lower bound by block-coordinate ascent with exact block LPs.

    Rows with positive multipliers are stored.  When a block is solved, the
    stored rows lying entirely inside it are released (their contribution is
    added back to r) and re-optimised together with newly separated rows, so a
    block step is an exact maximisation of (*) over all multipliers of rows
    inside the block."""

    def __init__(self, g: Graph, sup: Support | None = None):
        self.g = g
        self.sup = sup if sup is not None else build_support(g)
        s = self.sup
        self.ptr, self.idx, self.pid = _pgraph(g.n, g.indptr, g.indices, s.eid, s.n2ptr,
                                               s.n2idx, s.n2id)
        self.r = np.concatenate([np.ones(s.npos), -np.ones(s.nneg)])
        self.K = float(s.nneg)
        self.by = 0.0
        self.history = []
        self.t0 = time.time()
        # stored rows (global pair ids), grown in chunks
        self._chunks = []   # list of dicts with ptr, idx, val, b, y, anchor, alive
        self._anchor_index = None

    def bound(self) -> float:
        return self.K + self.by + float(np.minimum(self.r, 0.0).sum())

    def certified(self) -> int:
        return int(np.ceil(self.bound() - 1e-6))

    # -- stored rows -------------------------------------------------------
    def _rows_inside(self, pos):
        """Collect alive stored rows whose pairs all lie inside the block."""
        out = []
        pu, pv = self.sup.pu, self.sup.pv
        for ci, ch in enumerate(self._chunks):
            cand = np.flatnonzero(ch["alive"] & (pos[ch["anchor"]] >= 0))
            if len(cand) == 0:
                continue
            ptr, idx = ch["ptr"], ch["idx"]
            fl, seg = _gather(ptr, cand)
            flat = idx[fl]
            ok = (pos[pu[flat]] >= 0) & (pos[pv[flat]] >= 0)
            bad = np.zeros(len(cand), dtype=bool)
            np.logical_or.at(bad, seg, ~ok)
            sel = cand[~bad]
            if len(sel):
                out.append((ci, sel))
        return out

    def _store(self, rptr, gidx, rval, rb, y, anchor):
        keep = y > 1e-12
        if not keep.any():
            return
        sel = np.flatnonzero(keep)
        lens = rptr[sel + 1] - rptr[sel]
        nptr = np.zeros(len(sel) + 1, dtype=np.int64)
        nptr[1:] = np.cumsum(lens)
        flat, _ = _gather(rptr, sel)
        self._chunks.append({"ptr": nptr, "idx": gidx[flat], "val": rval[flat],
                             "b": rb[sel], "y": y[sel], "anchor": anchor[sel],
                             "alive": np.ones(len(sel), dtype=bool)})

    @property
    def nrows(self) -> int:
        return int(sum(ch["alive"].sum() for ch in self._chunks))

    # -- block step --------------------------------------------------------
    def solve_block(self, nodes: np.ndarray, reuse: bool = True, **kw) -> float:
        nodes = np.asarray(nodes, dtype=np.int64)
        pos = np.full(self.g.n, -1, dtype=np.int64)
        pos[nodes] = np.arange(len(nodes))
        lptr, lidx, lgid = _local_pgraph(nodes, pos, self.ptr, self.idx, self.pid)
        if lidx.shape[0] == 0:
            return 0.0
        gids, lpid = np.unique(lgid, return_inverse=True)
        old_bound = self.bound()
        released = self._rows_inside(pos) if reuse else []
        # release stored rows: r += A^T y, by -= b^T y
        r_save = self.r[gids].copy()
        by_save = self.by
        init = []
        for ci, sel in released:
            ch = self._chunks[ci]
            flat, rowof = _gather(ch["ptr"], sel)
            np.add.at(self.r, ch["idx"][flat], ch["val"][flat] * ch["y"][sel][rowof])
            self.by -= float(ch["b"][sel] @ ch["y"][sel])
            init.append((ch, sel, flat))
        costs = self.r[gids].copy()
        lp = _BlockLP(costs, len(nodes), lptr, lidx, lpid.astype(np.int64),
                      npos=int(np.searchsorted(gids, self.sup.npos)))
        # initial rows (translated to local pair ids)
        for ch, sel, flat in init:
            lens = ch["ptr"][sel + 1] - ch["ptr"][sel]
            starts = np.zeros(len(sel), dtype=np.int64)
            starts[1:] = np.cumsum(lens)[:-1]
            loc = np.searchsorted(gids, ch["idx"][flat])
            lp._push(starts, loc, ch["val"][flat], ch["b"][sel])
        xloc = lp.solve(**kw)
        self.last_block = {"gids": gids, "x": np.clip(xloc, 0.0, 1.0),
                           "converged": bool(lp.converged), "whole": len(nodes) == self.g.n}
        y = lp.duals()
        rptr, ridx, rval, rb = lp.rows()
        rloc = costs.copy()
        if len(rb):
            _apply_duals(rptr, ridx, rval, y, rloc)
        by = float(rb @ y) if len(rb) else 0.0
        # tentative bound
        self.r[gids] = rloc
        self.by += by
        if self.bound() + 1e-7 < old_bound:
            # revert (numerical trouble or time limit): restore stored rows
            self.r[gids] = r_save
            self.by = by_save
            return 0.0
        for ch, sel, _ in init:
            ch["alive"][sel] = False
        if len(rb):
            gidx = gids[ridx]
            # anchor = smallest endpoint of the first pair of the row
            first = gidx[rptr[:-1]]
            anchor = np.minimum(self.sup.pu[first], self.sup.pv[first]).astype(np.int64)
            self._store(rptr, gidx, rval, rb, y, anchor)
        self._compact()
        return self.bound() - old_bound

    def _compact(self):
        if len(self._chunks) > 64:
            chs = [c for c in self._chunks if c["alive"].any()]
            self._chunks = chs

    def sweep(self, labels: np.ndarray, **kw):
        """Solve every block of the partition ``labels`` once."""
        order = np.argsort(labels, kind="stable")
        ls = labels[order]
        cuts = np.flatnonzero(np.diff(ls)) + 1
        for blk in np.split(order, cuts):
            if len(blk) >= 3:
                self.solve_block(blk, **kw)
        self.history.append((time.time() - self.t0, self.bound()))
        return self.bound()


def metis_blocks(g: Graph, block_size: int, seed: int = 0) -> np.ndarray:
    """Vertex partition into ~n/block_size parts minimising the positive edge cut.

    Different seeds give different partitions (the vertices are randomly
    relabelled before calling METIS)."""
    import pymetis
    from .graph import from_edges
    k = max(1, int(np.ceil(g.n / block_size)))
    if k == 1:
        return np.zeros(g.n, dtype=np.int64)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(g.n)
    h = from_edges(g.n, perm[g.edges()])
    opts = pymetis.Options(seed=int(seed))
    adj = pymetis.CSRAdjacency(h.indptr, h.indices) if hasattr(pymetis, 'CSRAdjacency') else None
    if adj is not None:
        _, parts = pymetis.part_graph(k, adjacency=adj, options=opts)
    else:
        _, parts = pymetis.part_graph(k, xadj=h.indptr, adjncy=h.indices, options=opts)
    parts = np.asarray(parts, dtype=np.int64)
    out = np.empty(g.n, dtype=np.int64)
    out[:] = parts[perm]
    return out


def gap_blocks(bd: "BlockDualBound", labels: np.ndarray, size: int, time_limit: float,
               rng=None, block_kw=None, verbose=False, min_gain: float = 1e-3):
    """Dual-gap guided block steps: repeatedly grow a block of ``size`` vertices
    (BFS in G+, preferring high-gap vertices) around the vertex with the largest
    share of cost(C) - LB and re-solve it.  Returns the number of blocks."""
    from .lns import gap_map, grow_block
    rng = np.random.default_rng(rng)
    block_kw = block_kw or {}
    g = bd.g
    t0 = time.time()
    used = np.zeros(g.n, dtype=np.int64)
    it = 0
    while time.time() - t0 < time_limit:
        vg, _, _ = gap_map(bd, labels)
        eff = vg / (1.0 + used)
        seed = int(np.argmax(eff))
        if eff[seed] <= 1e-9:
            break
        B = grow_block(g, seed, size, vg, rng)
        used[B] += 1
        before = bd.bound()
        bd.solve_block(B, time_limit=max(1.0, min(block_kw.get("time_limit", 120),
                                                   time_limit - (time.time() - t0))),
                       **{k: v for k, v in block_kw.items() if k != "time_limit"})
        it += 1
        if verbose:
            print(f"  gap block {it}: {before:.1f} -> {bd.bound():.1f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    bd.history.append((time.time() - bd.t0, bd.bound()))
    return it


def add_star_packing(bd: "BlockDualBound", rp, rpairs, rk):
    """Install a star packing as initial multipliers (y = 1 on each star row).

    Row in >= form:  sum_t x_vt - sum_{tt'} x_tt' >= (k - 1) - k(k - 1)/2 ."""
    ns = len(rk)
    if ns == 0:
        return
    lens = rp[1:] - rp[:-1]
    k = rk.astype(np.float64)
    vals = np.empty(len(rpairs))
    # first k entries of each row are centre pairs (+1), the rest leaf pairs (-1)
    pos_in_row = np.arange(len(rpairs)) - np.repeat(rp[:-1], lens)
    vals[:] = np.where(pos_in_row < np.repeat(rk, lens), 1.0, -1.0)
    b = (k - 1.0) - k * (k - 1.0) / 2.0
    y = np.ones(ns)
    np.subtract.at(bd.r, rpairs, vals)
    bd.by += float(b.sum())
    first = rpairs[rp[:-1]]
    anchor = np.minimum(bd.sup.pu[first], bd.sup.pv[first]).astype(np.int64)
    bd._store(rp.astype(np.int64), rpairs.astype(np.int64), vals, b, y, anchor)
