"""Small-core test.  For a tight instance (OPT(G)=|H|) and the set A of anchorable vertices (all incident edges packed),
S subset of A is a 'face' if some optimal solution makes every vertex of S a singleton (=> all tunnel locks on S can be intact at cost 0).
2-decomposability test: if every pair of S is a face, is S a face?  (Minimal non-face of size >=3 would be a long lock-core.)"""
import sys, itertools, random, numpy as np, networkx as nx
from w6c_lib import *

def max_packing_ilp(A, rng=None, nsol=1):
    import scipy.sparse as sp
    from scipy.optimize import milp, LinearConstraint, Bounds
    P3 = p3_list(A)
    pairs = sorted({q for t in P3 for q in pairs_of(t)}); idx = {p: i for i, p in enumerate(pairs)}
    rows = []; cols = []
    for k, t in enumerate(P3):
        for q in pairs_of(t): rows.append(idx[q]); cols.append(k)
    M = sp.csr_matrix(([1] * len(rows), (rows, cols)), shape=(len(pairs), len(P3)))
    c = -np.ones(len(P3))
    if rng is not None: c = c - rng.random(len(P3)) * 0.1
    r = milp(c, constraints=LinearConstraint(M, -np.inf, 1), integrality=np.ones(len(P3)), bounds=Bounds(0, 1))
    return [P3[k] for k in range(len(P3)) if r.x[k] > 0.5]

def face(A, S):
    fixed = {}
    n = len(A)
    for a in S:
        for t in range(n):
            if t != a: fixed[(a, t) if a < t else (t, a)] = 0
    return ilp_opt(A, fixed=fixed)

def test(A, trials=5, seed=0, verbose=False):
    rng = np.random.default_rng(seed)
    n = len(A)
    opt, _ = ilp_opt(A)
    res = []
    for _ in range(trials):
        H = max_packing_ilp(A, rng)
        if len(H) != opt: return None
        used = set()
        for t in H: used |= pairs_of(t)
        anch = [v for v in range(n) if all(key(v, u) in used for u in range(n) if A[v][u])]
        # singleton feasible individually and pairwise
        ok1 = [a for a in anch if (face(A, [a]) or (None, 99))[0] == opt]
        comp = nx.Graph(); comp.add_nodes_from(ok1)
        for a, b in itertools.combinations(ok1, 2):
            r = face(A, [a, b])
            if r and r[0] == opt: comp.add_edge(a, b)
        worst = 0; viol = None
        for cl in nx.find_cliques(comp):
            if len(cl) < 3: continue
            r = face(A, cl)
            if not (r and r[0] == opt):
                viol = cl; break
        res.append((len(anch), len(ok1), viol))
        if viol: return (H, anch, viol)
    return res
