import sys, itertools, random, time
import numpy as np
sys.argv = ['x', '5']
from w11_maj5 import *
from pysat.solvers import Cadical153
rnd = random.Random(1); np.random.seed(1)
pats = pat_vars()
def solve_single(R, NO_=5):
    solver = Cadical153(); var = {}
    for p in pats:
        for o in range(NO_): var[(p, o)] = len(var) + 1
        solver.add_clause([var[(p, o)] for o in range(NO_)])
        for o1 in range(NO_):
            for o2 in range(o1 + 1, NO_): solver.add_clause([-var[(p, o1)], -var[(p, o2)]])
    def base_out(a, b, c):
        if a == b or a == c: return a
        if b == c: return b
        return None
    U, V, W = check_relation(R, None)
    k = R.shape[1]
    Rset = set(map(tuple, R.tolist()))
    # enumerate all combos of patterns (small k): constraint = forbid pattern-value combos whose output not in R
    for u, v, w in zip(U.tolist(), V.tolist(), W.tolist()):
        ps = [(u[c], v[c], w[c]) for c in range(k)]
        free = [c for c in range(k) if base_out(*ps[c]) is None]
        if not free: 
            out = tuple(base_out(*ps[c]) for c in range(k))
            if out not in Rset: return False, ('fixed', u, v, w)
            continue
        for vals in itertools.product(range(NO_), repeat=len(free)):
            out = [base_out(*ps[c]) for c in range(k)]
            for c, o in zip(free, vals): out[c] = o
            if tuple(out) not in Rset:
                solver.add_clause([-var[(ps[c], o)] for c, o in zip(free, vals)])
    return solver.solve(), None
bad = []
t0 = time.time()
cnt = 0
for name, cf in (('link', configs_link()),):
    seen = set()
    for (mem, nv, hard, X) in cf:
        R = rel5(mem, nv, hard, X)
        key = R.tobytes()
        if key in seen or len(R) == 0: continue
        seen.add(key)
        ok, why = solve_single(R)
        cnt += 1
        if not ok: bad.append((name, mem, hard, why))
print('link relations tested', cnt, 'individually UNSAT:', len(bad))
for b in bad[:3]: print(b[:3], b[3])

bad=[]; cnt=0
seen=set()
for (mem, nv, hard, X) in configs_k3():
    R = rel5(mem, nv, hard, X)
    key = R.tobytes()
    if key in seen or len(R) == 0: continue
    seen.add(key)
    ok, why = solve_single(R)
    cnt += 1
    if not ok: bad.append((mem, hard, why))
print('k3 relations tested', cnt, 'individually UNSAT:', len(bad), round(time.time()-t0))
for b in bad[:3]: print(b[0], b[1], b[2])
