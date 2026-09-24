"""Independent check of a lower-bound certificate.

A certificate is a list of rows  sum_p a_rp x_p >= b_r  over vertex pairs p
with multipliers y_r >= 0 (written by BlockDualBound.certificate).  The check
uses only the graph and the certificate, not the solver code:

1. every pair of a row is at distance at most two (an edge or a pair with a
   common neighbour), so it belongs to the relaxation P;
2. every row is valid for every clustering that separates the far pairs:
   rows on at most MAX_BRUTE vertices by enumerating all such partitions of
   their vertex set, larger rows only if they are star rows, for which the
   minimum has a closed form;
3. y is rounded down to multiples of 2^-SCALE (this keeps y >= 0, so the rows
   stay a valid certificate), and

       LB = b^T y + sum_p [min(0, c_p - (A^T y)_p) - min(0, c_p)]

   is evaluated in exact integer arithmetic (c_p = +1 on edges, -1 otherwise).

For every clustering C that separates far pairs, cost(C) = |N2| + sum_P c_p x_p
>= LB by weak duality; optimal clusterings separate far pairs, so LB <= OPT
and ceil(LB) is a certified lower bound.

usage: check_certificate.py GRAPH CERT.npz [GRAPH CERT.npz ...]
"""
import sys
from fractions import Fraction
from math import ceil

import numba as nb
import numpy as np

SCALE = 30
MAX_BRUTE = 9


@nb.njit(cache=True)
def _adjacent(indptr, indices, u, v):
    lo, hi = indptr[u], indptr[u + 1]
    while lo < hi:
        mid = (lo + hi) // 2
        if indices[mid] < v:
            lo = mid + 1
        else:
            hi = mid
    return lo < indptr[u + 1] and indices[lo] == v


@nb.njit(cache=True)
def _close(indptr, indices, u, v):
    """Adjacent or with a common neighbour (neighbour lists are sorted)."""
    if _adjacent(indptr, indices, u, v):
        return True
    i, j = indptr[u], indptr[v]
    while i < indptr[u + 1] and j < indptr[v + 1]:
        a, b = indices[i], indices[j]
        if a == b:
            return True
        if a < b:
            i += 1
        else:
            j += 1
    return False


@nb.njit(cache=True)
def _check_pairs(indptr, indices, u, v, out_sign):
    """Sign of every pair (+1 edge, -1 distance two); returns the number of bad pairs."""
    bad = 0
    for k in range(u.shape[0]):
        a, b = u[k], v[k]
        if a == b:
            bad += 1
            continue
        if _adjacent(indptr, indices, a, b):
            out_sign[k] = 1
        elif _close(indptr, indices, a, b):
            out_sign[k] = -1
        else:
            bad += 1
    return bad


@nb.njit(cache=True)
def _brute_rows(indptr, indices, ptr, u, v, val, b, rows, max_brute):
    """For each row in ``rows`` (at most max_brute vertices): minimum of the
    left-hand side over all partitions of its vertex set that separate its far
    pairs; returns a flag per row (1 valid, 0 invalid, -1 too large)."""
    ok = np.zeros(rows.shape[0], dtype=np.int64)
    U = np.empty(max_brute, dtype=np.int64)
    li = np.empty(0, dtype=np.int64)
    far = np.zeros((max_brute, max_brute), dtype=np.bool_)
    rgs = np.zeros(max_brute, dtype=np.int64)
    mx = np.zeros(max_brute, dtype=np.int64)
    for q in range(rows.shape[0]):
        r = rows[q]
        s, e = ptr[r], ptr[r + 1]
        k = 0
        too_big = False
        for p in range(s, e):
            for w in (u[p], v[p]):
                found = False
                for t in range(k):
                    if U[t] == w:
                        found = True
                        break
                if not found:
                    if k == max_brute:
                        too_big = True
                        break
                    U[k] = w
                    k += 1
            if too_big:
                break
        if too_big:
            ok[q] = -1
            continue
        li = np.empty(2 * (e - s), dtype=np.int64)
        for p in range(s, e):
            for t in range(k):
                if U[t] == u[p]:
                    li[2 * (p - s)] = t
                if U[t] == v[p]:
                    li[2 * (p - s) + 1] = t
        for a in range(k):
            for c in range(a + 1, k):
                f = not _close(indptr, indices, min(U[a], U[c]), max(U[a], U[c]))
                far[a, c] = f
                far[c, a] = f
        # enumerate restricted growth strings
        best = 1e18
        for t in range(k):
            rgs[t] = 0
            mx[t] = 0
        while True:
            feasible = True
            for a in range(k):
                for c in range(a + 1, k):
                    if far[a, c] and rgs[a] == rgs[c]:
                        feasible = False
                        break
                if not feasible:
                    break
            if feasible:
                lhs = 0.0
                for p in range(s, e):
                    if rgs[li[2 * (p - s)]] != rgs[li[2 * (p - s) + 1]]:
                        lhs += val[p]
                if lhs < best:
                    best = lhs
            # next restricted growth string
            t = k - 1
            while t > 0 and rgs[t] > mx[t - 1]:
                t -= 1
            if t == 0:
                break
            rgs[t] += 1
            m = max(mx[t - 1], rgs[t])
            mx[t] = m
            for z in range(t + 1, k):
                rgs[z] = 0
                mx[z] = m
        ok[q] = 1 if best >= b[r] - 1e-9 else 0
    return ok


def _max_independent(adj):
    """Exact maximum independent set size of a small graph (dict of sets)."""
    best = 0

    def rec(cand, size):
        nonlocal best
        if size + len(cand) <= best:
            return
        if not cand:
            best = max(best, size)
            return
        v = max(cand, key=lambda x: len(adj[x] & cand))
        rec(cand - {v} - adj[v], size + 1)
        if adj[v] & cand:
            rec(cand - {v}, size)

    rec(set(adj), 0)
    return best


def _star_valid(g, rows_u, rows_v, rows_val, b, close):
    """Row  sum_t x_vt - sum_{tt' in R} x_tt' >= b  with centre v and leaves T.
    Minimum over far-separating clusterings: |T| - |R| - M with
    M = max over S (S + v far-free) of |S| - e_R(S).  If R contains every
    non-far pair of T, a far-free S is a clique of R and M <= 1; otherwise
    M <= alpha(R[T]) (a component of R[S] with s vertices has >= s - 1 edges)."""
    plus = [(a, c) for a, c, x in zip(rows_u, rows_v, rows_val) if x == 1.0]
    minus = [(a, c) for a, c, x in zip(rows_u, rows_v, rows_val) if x == -1.0]
    if len(plus) + len(minus) != len(rows_val) or not plus:
        return False
    common = set(plus[0])
    for p in plus:
        common &= set(p)
    if len(common) != 1:
        return False
    v = common.pop()
    T = [a if c == v else c for a, c in plus]
    if len(set(T)) != len(T):
        return False
    Ts = set(T)
    R = set()
    for a, c in minus:
        if a not in Ts or c not in Ts or a == c:
            return False
        R.add((min(a, c), max(a, c)))
    if len(R) != len(minus):
        return False
    Tl = sorted(Ts)
    complete = all((a, c) in R or not close(a, c)
                   for i, a in enumerate(Tl) for c in Tl[i + 1:])
    if complete:
        M = 1
    else:
        adj = {t: set() for t in Tl}
        for a, c in R:
            adj[a].add(c)
            adj[c].add(a)
        M = _max_independent(adj)
    return b <= len(T) - len(R) - M + 1e-9


def check(g, cert):
    ptr, u, v = cert["ptr"].astype(np.int64), cert["u"].astype(np.int64), cert["v"].astype(np.int64)
    val, b, y = cert["val"].astype(np.float64), cert["b"].astype(np.float64), cert["y"].astype(np.float64)
    nrows = len(b)
    report = {"rows": nrows, "nnz": len(u)}
    lo, hi = np.minimum(u, v), np.maximum(u, v)
    if (y < 0).any() or not np.all(np.isfinite(y)):
        raise ValueError("negative or non-finite multiplier")
    if not (np.all(val == np.round(val)) and np.all(b == np.round(b))):
        raise ValueError("non-integral row data")
    sign = np.zeros(len(u), dtype=np.int64)
    bad = _check_pairs(g.indptr, g.indices, lo, hi, sign)
    if bad:
        raise ValueError(f"{bad} row entries are not pairs of P")
    # validity of every row
    flags = _brute_rows(g.indptr, g.indices, ptr, lo, hi, val, b,
                        np.arange(nrows, dtype=np.int64), MAX_BRUTE)
    if (flags == 0).any():
        raise ValueError(f"{int((flags == 0).sum())} rows are not valid")
    big = np.flatnonzero(flags == -1)

    def close(a, c):
        return bool(_close(g.indptr, g.indices, min(a, c), max(a, c)))

    for r in big:
        s, e = ptr[r], ptr[r + 1]
        if not _star_valid(g, lo[s:e].tolist(), hi[s:e].tolist(), val[s:e].tolist(), b[r], close):
            raise ValueError(f"row {r} ({e - s} entries) could not be verified")
    report["brute_rows"] = int((flags == 1).sum())
    report["star_rows"] = len(big)
    # exact evaluation
    yi = np.floor(y * 2.0 ** SCALE).astype(np.int64)
    key = lo * g.n + hi
    uniq, inv = np.unique(key, return_inverse=True)
    rowof = np.repeat(np.arange(nrows), np.diff(ptr))
    a = np.zeros(len(uniq), dtype=np.int64)
    np.add.at(a, inv, val.astype(np.int64) * yi[rowof])
    c = np.zeros(len(uniq), dtype=np.int64)
    c[inv] = sign
    one = 1 << SCALE
    terms = np.minimum(0, c * one - a) - np.minimum(0, c * one)
    total = sum(int(x) for x in b.astype(np.int64).astype(object) * yi.astype(object))
    total += sum(int(x) for x in terms.astype(object))
    lb = Fraction(total, one)
    report["lb_exact"] = float(lb)
    report["certified"] = ceil(lb)
    report["solver_bound"] = float(cert["bound"])
    return report


def main(argv):
    sys.path.insert(0, ".")
    sys.path.insert(0, "experiments")
    import datasets as D
    for name, path in zip(argv[0::2], argv[1::2]):
        g = D.load(name)
        rep = check(g, np.load(path))
        print(name, rep, flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
