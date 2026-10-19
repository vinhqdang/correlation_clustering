"""Wave 5 B: ILP spot checks (n = 11..13) of EX1/EX2 (endpoint hard degree, exceptional cluster size).
For each instance (pair-maximal H, budget l): for each vertex u,
  (EX1) if d_h(u) > D + l then NO solution of cost <= |H|+l edits a hard pair at u   [ILP infeasible]
  (EX2) max |cluster(u)| over solutions of cost <= |H|+l editing a hard pair at u  <=  2D+2l+1.
usage: python3 -I w5b_ilp.py N SEED NMIN NMAX"""
import sys, random, itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix
sys.path.insert(0, '.')
from w5b_lib import *
_TR = {}
def trans(n):
    if n not in _TR:
        P = {p: i for i, p in enumerate(prs(n))}; r = []; c = []; v = []; k = 0
        for a, b, d in itertools.combinations(range(n), 3):
            ab, bd, ad = P[(a, b)], P[(b, d)], P[(a, d)]
            for x, y, z in ((ab, bd, ad), (ab, ad, bd), (bd, ad, ab)):
                r += [k, k, k]; c += [x, y, z]; v += [1, 1, -1]; k += 1
        _TR[n] = coo_matrix((v, (r, c)), shape=(k, len(P))).tocsr()
    return _TR[n]
def solve(A, budget, hard_at, objpairs, sense):
    n = len(A); ps = prs(n); m = len(ps); g = np.array([A[a][b] for a, b in ps]); s = 1 - 2 * g
    T = trans(n); M = T.multiply(s[None, :]).tocsr()
    cons = [LinearConstraint(M, -np.inf, 1 - T @ g), LinearConstraint(np.ones((1, m)), 0, budget)]
    row = np.zeros(m)
    for i, (a, b) in enumerate(ps):
        if hard_at(a, b): row[i] = 1
    cons.append(LinearConstraint(row.reshape(1, -1), 1, np.inf))
    obj = np.zeros(m)
    for i, (a, b) in enumerate(ps):
        if objpairs(a, b): obj[i] = s[i] * sense   # y = g + s e ; maximise sum y -> maximise s*e
    r = milp(-obj, constraints=cons, integrality=np.ones(m), bounds=Bounds(0, 1))
    if r.status == 2: return None
    if r.status != 0: raise RuntimeError(r.message)
    e = np.round(r.x); return float(sum(obj[i] * e[i] for i in range(m)))
def run(ninst, seed, nmin, nmax):
    rng = random.Random(seed); cnt = dict(inst=0, ex1_checked=0, ex1_feasible_wrongly=0, ex2_checked=0, ex2_viol=0, max_slack=-99, ex1_nontrivial=0)
    for it in range(ninst):
        n = rng.randint(nmin, nmax); A = rand_graph(n, rng)
        if A.sum() == 0: continue
        H = pair_maximal_packing(A, rng)
        if not H: continue
        pk = soft_matrix(n, H); D = int(pk.sum(1).max()); comp, hg = hard_comps(A, pk); dh = hg.sum(1)
        # OPT via ILP with no 'hard_at' requirement
        ps = prs(n); m = len(ps); g = np.array([A[a][b] for a, b in ps]); s = 1 - 2 * g; T = trans(n)
        r = milp(np.ones(m), constraints=[LinearConstraint(T.multiply(s[None, :]).tocsr(), -np.inf, 1 - T @ g)],
                 integrality=np.ones(m), bounds=Bounds(0, 1)); opt = int(round(r.fun)); l = opt - len(H)
        cnt['inst'] += 1
        for ll in sorted(set([l, l + 1])):
            budget = len(H) + ll
            for u in range(n):
                hard_at = lambda a, b, u=u: (a == u or b == u) and not pk[a, b]
                if dh[u] == 0 and not any(not pk[u, w] for w in range(n) if w != u): continue
                if dh[u] > D + ll:
                    cnt['ex1_checked'] += 1; cnt['ex1_nontrivial'] += 1
                    if solve(A, budget, hard_at, lambda a, b: False, 1) is not None: cnt['ex1_feasible_wrongly'] += 1
                else:
                    # EX2: maximise sum_w y_uw  (cluster size - 1)
                    val = solve(A, budget, hard_at, lambda a, b, u=u: a == u or b == u, 1)
                    if val is None: continue
                    size = 1 + int(round(val + sum(g[i] for i, (a, b) in enumerate(ps) if a == u or b == u)))
                    cnt['ex2_checked'] += 1; cnt['max_slack'] = max(cnt['max_slack'], size - (2 * D + 2 * ll + 1))
                    if size > 2 * D + 2 * ll + 1: cnt['ex2_viol'] += 1
    return cnt
if __name__ == '__main__':
    print(run(*map(int, sys.argv[1:5])))
