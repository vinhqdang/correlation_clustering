"""Wave 8 library (written fresh; imports only numpy/pysat/itertools).

Conventions.  Vertices 0..n-1.  A member is a triple (a,b,c) with centre b: ab, bc in E, ac not in E.
A solution is a partition; ZE (zero excess) = exactly one edited pair per member and no other edited pair.
Togetherness variable t(u,v).  All SAT work is done with pysat (CaDiCaL).
"""
import itertools
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType


def key(a, b):
    return (a, b) if a < b else (b, a)


def pairs_of(t):
    a, b, c = t
    return {key(a, b), key(b, c), key(a, c)}


class VarPool:
    def __init__(self):
        self.n = 0

    def new(self):
        self.n += 1
        return self.n


def ze_cnf(n, E, members, vp, singletons=(), fixed_apart=(), fixed_together=()):
    """Return (clauses, T) where T[(u,v)] is the variable for 'u,v together' (u<v).
    E: set of edges (u,v) u<v.  members: list of triples.  ZE conditions are encoded."""
    T = {}
    for u, v in itertools.combinations(range(n), 2):
        T[(u, v)] = vp.new()
    cl = []
    for u, v, w in itertools.combinations(range(n), 3):
        a, b, c = T[(u, v)], T[(v, w)], T[(u, w)]
        cl += [[-a, -b, c], [-a, -c, b], [-b, -c, a]]
    soft = set()
    for m in members:
        soft |= pairs_of(m)
    for p in T:
        if p not in soft:
            cl.append([T[p]] if p in E else [-T[p]])
    for (a, b, c) in members:
        # edited(ab) = not t(ab) ; edited(bc) = not t(bc) ; edited(ac) = t(ac)
        e1, e2, e3 = -T[key(a, b)], -T[key(b, c)], T[key(a, c)]
        cl += [[e1, e2, e3], [-e1, -e2], [-e1, -e3], [-e2, -e3]]
    for x in singletons:
        for v in range(n):
            if v != x:
                cl.append([-T[key(x, v)]])
    for (u, v) in fixed_apart:
        cl.append([-T[key(u, v)]])
    for (u, v) in fixed_together:
        cl.append([T[key(u, v)]])
    return cl, T


def ze_solve(n, E, members, singletons=(), fixed_apart=(), fixed_together=(), all_solutions=False, limit=None):
    """Return a ZE partition (labels list) with the given singletons, or None.  With all_solutions
    return the list of all (up to limit)."""
    vp = VarPool()
    cl, T = ze_cnf(n, E, members, vp, singletons, fixed_apart, fixed_together)
    s = Cadical153(bootstrap_with=cl)
    out = []
    while s.solve():
        model = set(l for l in s.get_model() if l > 0)
        lab = [-1] * n
        k = 0
        for v in range(n):
            if lab[v] < 0:
                lab[v] = k
                for w in range(v + 1, n):
                    if T[(v, w)] in model:
                        lab[w] = k
                k += 1
        if not all_solutions:
            s.delete()
            return lab
        out.append(lab)
        if limit and len(out) >= limit:
            break
        s.add_clause([(-T[p] if T[p] in model else T[p]) for p in T])
    s.delete()
    return out if all_solutions else None


def partition_cost(n, E, lab):
    return sum(1 for u, v in itertools.combinations(range(n), 2) if ((u, v) in E) != (lab[u] == lab[v]))


def check_instance(n, E, members):
    """validity: members induced P3s, modification-disjoint."""
    seen = set()
    for (a, b, c) in members:
        if not (key(a, b) in E and key(b, c) in E and key(a, c) not in E):
            return False, ('notP3', (a, b, c))
        ps = pairs_of((a, b, c))
        if ps & seen:
            return False, ('overlap', (a, b, c))
        seen |= ps
    return True, None


def ilp_opt_ok(n, E, members, extra_vertices=0):
    """exact OPT by ILP (HiGHS) for cross-checking small instances; returns OPT."""
    import numpy as np
    from scipy.optimize import milp, LinearConstraint, Bounds
    import scipy.sparse as sp
    P = list(itertools.combinations(range(n), 2))
    idx = {p: i for i, p in enumerate(P)}
    m = len(P)
    e = np.array([1.0 if p in E else 0.0 for p in P])
    c = 1 - 2 * e
    r, cc, vv = [], [], []
    k = 0
    for i, j, l in itertools.combinations(range(n), 3):
        a, b, d = idx[(i, j)], idx[(j, l)], idx[(i, l)]
        for (p, q, s_) in ((a, b, d), (a, d, b), (b, d, a)):
            r += [k, k, k]
            cc += [p, q, s_]
            vv += [1, 1, -1]
            k += 1
    A = sp.csr_matrix((vv, (r, cc)), shape=(k, m))
    res = milp(c, constraints=[LinearConstraint(A, -np.inf, 1)], integrality=np.ones(m), bounds=Bounds(0, 1))
    return int(round(e.sum() + res.fun))


def labels_to_clusters(lab):
    d = {}
    for v, l in enumerate(lab):
        d.setdefault(l, []).append(v)
    return sorted(d.values())
