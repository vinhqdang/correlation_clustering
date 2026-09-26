"""Independent check of a lower-bound certificate.

A certificate is a list of rows  sum_p a_rp x_p >= b_r  over vertex pairs p
with multipliers y_r >= 0 (written by BlockDualBound.certificate).  The check
uses only the raw instance file and the certificate.  It imports nothing from
the solver package (ccbench) or its data loader: the instance is parsed here,
from the file as distributed by SNAP or PACE.

0. the raw file is parsed by the reader below (SNAP edge list, optionally
   signed, or PACE .gr); its SHA-256, n, m and a hash of the canonical edge
   list are compared with the values recorded in the certificate, so a
   certificate cannot silently refer to a different graph or a different
   vertex numbering;
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

Instance semantics (Section 2 of the paper): the vertices of a SNAP instance
are the node ids that occur in a kept edge, numbered 0..n-1 in increasing id
order; in a signed file only rows with a positive value in the sign column are
kept; edges are undirected, self-loops and duplicates are dropped.  PACE .gr
files are 1-indexed with a ``p cep n m`` header.

The upper bound is re-checked in the same way: ``clustering_cost`` evaluates
the cost of an archived clustering (one label per vertex, in the numbering
above) on the instance parsed here.

usage: check_certificate.py RAW_FILE CERT.npz [RAW_FILE CERT.npz ...]
       check_certificate.py --dir DATA_DIR CERT_DIR OUT.csv
       check_certificate.py --labels RAW_FILE FORMAT LABELS.npz
"""
import gzip
import hashlib
import sys
from fractions import Fraction
from math import ceil

import numba as nb
import numpy as np

SCALE = 30
MAX_BRUTE = 9


# --------------------------------------------------------------------------
# instance reader (deliberately independent of ccbench.graph)

class Instance:
    def __init__(self, n, u, v, raw_sha256):
        self.n = n
        self.m = len(u)
        self.raw_sha256 = raw_sha256
        # CSR with sorted neighbour lists
        src = np.concatenate([u, v])
        dst = np.concatenate([v, u])
        order = np.lexsort((dst, src))
        src, dst = src[order], dst[order]
        self.indptr = np.zeros(n + 1, dtype=np.int64)
        self.indptr[1:] = np.cumsum(np.bincount(src, minlength=n))
        self.indices = dst.astype(np.int64)
        self.edge_sha256 = edge_hash(n, u, v)


def edge_hash(n, u, v):
    """SHA-256 of n and the lexicographically sorted edge list (u < v, int64)."""
    order = np.lexsort((v, u))
    h = hashlib.sha256()
    h.update(np.int64(n).tobytes())
    h.update(np.ascontiguousarray(np.stack([u[order], v[order]], 1), dtype="<i8").tobytes())
    return h.hexdigest()


def _canonical(n, a, b):
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    keep = lo != hi
    key = np.unique(lo[keep] * n + hi[keep])
    return key // n, key % n


def read_instance(path, fmt, sign_col=-1):
    raw = open(path, "rb").read()
    sha = hashlib.sha256(raw).hexdigest()
    text = gzip.decompress(raw).decode() if path.endswith(".gz") else raw.decode()
    if fmt == "pace":
        n = None
        a, b = [], []
        for line in text.splitlines():
            if not line.strip() or line.startswith("c"):
                continue
            parts = line.split()
            if parts[0] == "p":
                n = int(parts[2])
                continue
            a.append(int(parts[0]) - 1)
            b.append(int(parts[1]) - 1)
        if n is None:
            raise ValueError("missing 'p cep' header")
        a, b = np.array(a, dtype=np.int64), np.array(b, dtype=np.int64)
        if len(a) and (min(a.min(), b.min()) < 0 or max(a.max(), b.max()) >= n):
            raise ValueError("vertex id out of range")
        u, v = _canonical(n, a, b)
        return Instance(n, u, v, sha)
    if fmt != "snap":
        raise ValueError(f"unknown format {fmt}")
    a, b = [], []
    for line in text.splitlines():
        if not line.strip() or line[0] in "#%":
            continue
        parts = line.replace(",", " ").split()
        if sign_col >= 0 and float(parts[sign_col]) <= 0:
            continue
        a.append(int(parts[0]))
        b.append(int(parts[1]))
    a, b = np.array(a, dtype=np.int64), np.array(b, dtype=np.int64)
    ids = np.unique(np.concatenate([a, b]))
    n = len(ids)
    u, v = _canonical(n, np.searchsorted(ids, a), np.searchsorted(ids, b))
    return Instance(n, u, v, sha)


# --------------------------------------------------------------------------
# pair and row checks

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
    pairs; returns a flag per row (1 valid, 0 invalid, -1 too large).  Row data
    are integers, so the sums are exact."""
    ok = np.zeros(rows.shape[0], dtype=np.int64)
    U = np.empty(max_brute, dtype=np.int64)
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
        best = np.int64(1) << 60
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
                lhs = np.int64(0)
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
        ok[q] = 1 if best >= b[r] else 0
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


def _star_valid(rows_u, rows_v, rows_val, b, close):
    """Row  sum_t x_vt - sum_{tt' in R} x_tt' >= b  with centre v and leaves T.
    Minimum over far-separating clusterings: |T| - |R| - M with
    M = max over S (S + v far-free) of |S| - e_R(S).  If R contains every
    non-far pair of T, a far-free S is a clique of R and M <= 1; otherwise
    M <= alpha(R[T]) (a component of R[S] with s vertices has >= s - 1 edges)."""
    plus = [(a, c) for a, c, x in zip(rows_u, rows_v, rows_val) if x == 1]
    minus = [(a, c) for a, c, x in zip(rows_u, rows_v, rows_val) if x == -1]
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
    return b <= len(T) - len(R) - M


def _meta(cert, key, default=None):
    if key not in cert:
        return default
    x = cert[key]
    return x.item() if hasattr(x, "item") else x


def check(g, cert):
    """g: Instance read by read_instance; cert: mapping of arrays."""
    report = {"n": g.n, "m": g.m, "raw_sha256": g.raw_sha256, "edge_sha256": g.edge_sha256}
    # 0. identity of the instance
    if int(_meta(cert, "n", g.n)) != g.n:
        raise ValueError(f"certificate is for n={_meta(cert, 'n')}, file has n={g.n}")
    for key, have in (("instance_m", g.m), ("edge_sha256", g.edge_sha256),
                      ("raw_sha256", g.raw_sha256)):
        want = _meta(cert, key)
        if want is not None and str(want) != str(have):
            raise ValueError(f"{key} mismatch: certificate {want}, file {have}")
    report["identity"] = "hash" if _meta(cert, "edge_sha256") is not None else "n only"
    ptr = np.asarray(cert["ptr"]).astype(np.int64)
    u, v = np.asarray(cert["u"]).astype(np.int64), np.asarray(cert["v"]).astype(np.int64)
    fval, fb = np.asarray(cert["val"], dtype=np.float64), np.asarray(cert["b"], dtype=np.float64)
    y = np.asarray(cert["y"], dtype=np.float64)
    nrows = len(fb)
    report.update({"rows": nrows, "nnz": len(u)})
    if len(ptr) != nrows + 1 or ptr[0] != 0 or ptr[-1] != len(u) or (np.diff(ptr) < 0).any():
        raise ValueError("malformed row pointer")
    if len(y) != nrows or (y < 0).any() or not np.all(np.isfinite(y)):
        raise ValueError("negative or non-finite multiplier")
    if not (np.all(fval == np.round(fval)) and np.all(fb == np.round(fb))):
        raise ValueError("non-integral row data")
    if len(u) and (min(u.min(), v.min()) < 0 or max(u.max(), v.max()) >= g.n):
        raise ValueError("vertex id out of range")
    val, b = fval.astype(np.int64), fb.astype(np.int64)
    lo, hi = np.minimum(u, v), np.maximum(u, v)
    # 1. pairs of P
    sign = np.zeros(len(u), dtype=np.int64)
    bad = _check_pairs(g.indptr, g.indices, lo, hi, sign)
    if bad:
        raise ValueError(f"{bad} row entries are not pairs of P")
    # 2. validity of every row
    flags = _brute_rows(g.indptr, g.indices, ptr, lo, hi, val, b,
                        np.arange(nrows, dtype=np.int64), MAX_BRUTE)
    if (flags == 0).any():
        raise ValueError(f"{int((flags == 0).sum())} rows are not valid")
    big = np.flatnonzero(flags == -1)

    def close(a, c):
        return bool(_close(g.indptr, g.indices, min(a, c), max(a, c)))

    for r in big:
        s, e = ptr[r], ptr[r + 1]
        if not _star_valid(lo[s:e].tolist(), hi[s:e].tolist(), val[s:e].tolist(), int(b[r]),
                           close):
            raise ValueError(f"row {r} ({e - s} entries) could not be verified")
    report["brute_rows"] = int((flags == 1).sum())
    report["star_rows"] = len(big)
    # 3. exact evaluation, in Python integers
    yi = np.floor(y * 2.0 ** SCALE).astype(np.int64)
    key = lo * g.n + hi
    uniq, inv = np.unique(key, return_inverse=True)
    rowof = np.repeat(np.arange(nrows), np.diff(ptr))
    # int64 accumulation is exact if no partial sum can exceed 2^62; this is
    # checked with an upper bound in floating point (with a generous margin)
    absbound = np.zeros(len(uniq))
    np.add.at(absbound, inv, np.abs(val) * yi[rowof].astype(np.float64))
    one = 1 << SCALE
    if len(absbound) and absbound.max() * 1.01 + one >= 2.0 ** 62:
        raise ValueError("coefficients too large for exact int64 accumulation")
    a = np.zeros(len(uniq), dtype=np.int64)
    np.add.at(a, inv, val * yi[rowof])
    c = np.zeros(len(uniq), dtype=np.int64)
    c[inv] = sign
    terms = np.minimum(0, c * one - a) - np.minimum(0, c * one)
    total = sum(int(bi) * int(yj) for bi, yj in zip(b.tolist(), yi.tolist()))
    total += sum(terms.tolist())
    lb = Fraction(total, one)
    report["lb_exact"] = float(lb)
    report["certified"] = ceil(lb)
    report["solver_bound"] = float(_meta(cert, "bound", float("nan")))
    return report


def clustering_cost(g, labels):
    """Exact number of disagreements of a clustering of the instance g."""
    labels = np.asarray(labels).astype(np.int64)
    if labels.shape != (g.n,):
        raise ValueError(f"{labels.shape[0]} labels for n={g.n}")
    src = np.repeat(np.arange(g.n), np.diff(g.indptr))
    inside = int((labels[src] == labels[g.indices]).sum()) // 2
    sizes = np.unique(labels, return_counts=True)[1]
    pairs_inside = sum(int(s) * (int(s) - 1) // 2 for s in sizes.tolist())
    return (g.m - inside) + (pairs_inside - inside)


# --------------------------------------------------------------------------

def check_file(raw, cert_path, fmt=None, sign_col=None):
    cert = dict(np.load(cert_path, allow_pickle=False))
    fmt = fmt or str(_meta(cert, "instance_format", "pace" if raw.endswith(".gr") else "snap"))
    sign_col = int(_meta(cert, "sign_col", -1)) if sign_col is None else sign_col
    return check(read_instance(raw, fmt, sign_col), cert)


def check_dir(data_dir, cert_dir, out_csv):
    """Check every certificate in cert_dir against the raw file named in it
    (key ``instance_file``, relative to data_dir) and write one CSV row each."""
    import csv
    import glob
    import os
    rows = []
    for path in sorted(glob.glob(os.path.join(cert_dir, "**", "*.npz"), recursive=True)):
        cert = dict(np.load(path, allow_pickle=False))
        if path.endswith(".labels.npz"):
            # a clustering without certificate (primal runs): check its cost only
            cert_path = path[:-len(".labels.npz")] + ".npz"
            if os.path.exists(cert_path):
                continue
            lab = dict(np.load(path, allow_pickle=False))
            row = {"certificate": os.path.relpath(path, cert_dir)}
            try:
                f = _meta(lab, "instance_file")
                if f is None:
                    raise ValueError("clustering names no instance file")
                g = read_instance(os.path.join(data_dir, f), str(_meta(lab, "instance_format")),
                                  int(_meta(lab, "sign_col", -1)))
                if str(_meta(lab, "edge_sha256")) != g.edge_sha256:
                    raise ValueError("edge_sha256 mismatch")
                row.update({"instance": f, "n": g.n, "m": g.m, "identity": "hash",
                            "cost": clustering_cost(g, lab["labels"]), "status": "ok",
                            "raw_sha256": g.raw_sha256, "edge_sha256": g.edge_sha256})
            except (ValueError, OSError) as exc:
                row["status"] = f"rejected: {exc}"
            rows.append(row)
            continue
        row = {"certificate": os.path.relpath(path, cert_dir)}
        try:
            f = _meta(cert, "instance_file")
            if f is None:
                raise ValueError("certificate names no instance file")
            row["instance"] = f
            row.update(check_file(os.path.join(data_dir, f), path))
            lab = path[:-4] + ".labels.npz"
            if os.path.exists(lab):
                g = read_instance(os.path.join(data_dir, f),
                                  str(_meta(cert, "instance_format")),
                                  int(_meta(cert, "sign_col", -1)))
                row["cost"] = clustering_cost(g, np.load(lab)["labels"])
            row["status"] = "ok"
        except (ValueError, OSError) as exc:
            row["status"] = f"rejected: {exc}"
        rows.append(row)
        print({k: row.get(k) for k in ("certificate", "status", "certified")}, flush=True)
    keys = ["certificate", "instance", "status", "identity", "n", "m", "rows", "nnz",
            "brute_rows", "star_rows", "lb_exact", "certified", "solver_bound", "cost", "raw_sha256",
            "edge_sha256"]
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    bad = [r for r in rows if r["status"] != "ok"]
    print(f"{len(rows) - len(bad)} of {len(rows)} certificates accepted", flush=True)
    return not bad


def main(argv):
    if argv and argv[0] == "--labels":
        g = read_instance(argv[1], argv[2])
        print(argv[3], clustering_cost(g, np.load(argv[3])["labels"]), flush=True)
        return
    if argv and argv[0] == "--dir":
        sys.exit(0 if check_dir(argv[1], argv[2], argv[3]) else 1)
    for raw, path in zip(argv[0::2], argv[1::2]):
        print(raw, check_file(raw, path), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
