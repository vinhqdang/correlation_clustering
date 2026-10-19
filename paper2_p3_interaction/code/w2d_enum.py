"""Exp 2b: enumerate ALL optimal edit sets (cap 400) for random instances, n=9..10, rand maximal H with l>0, |R|>=1.
Reports fraction of optima having an unpacked edit touching R / inside R (R-R). Seeds fixed."""
import random, itertools, json, collections
from w2d_common import *
def enum_opt(n, E, O, cap=400):
    ps = list(itertools.combinations(range(n), 2)); m = len(ps); cuts = []; sols = []
    while len(sols) < cap:
        r = solve(n, E, cap=O, extra=cuts)
        if r is None: return sols, True
        S = r[1]; sols.append(S)
        row = np.array([-1.0 if p in S else 1.0 for p in ps]); cuts.append((row, 1 - len(S), 1e9))
    return sols, False
rng = random.Random(20261019); tot = collections.Counter(); inst = 0; rows = []
for fam in ("gnp", "planted", "nearclique", "pathcliques", "blowup"):
    for n in (9, 10):
        for s in range(12):
            r2 = random.Random(f"enum-{fam}-{n}-{s}"); E = FAMILIES[fam](r2, n); O = opt(n, E); H = rand_maximal_packing(r2, n, E)
            l = O - len(H); V = vh(H)
            if l <= 0 or len(V) == n: continue
            PP = packed_pairs(H); sols, complete = enum_opt(n, E, O)
            if not complete: continue
            rt = sum(any(p not in PP and (p[0] not in V or p[1] not in V) for p in S) for S in sols)
            rr = sum(any(p not in PP and p[0] not in V and p[1] not in V for p in S) for S in sols)
            rows.append((fam, n, s, l, len(sols), rt, rr)); tot["inst"] += 1; tot["sols"] += len(sols); tot["rt"] += rt; tot["rr"] += rr
            tot["inst_all_rt"] += (rt == len(sols)); tot["inst_some_rt"] += (rt > 0); tot["inst_no_rt"] += (rt == 0)
print(dict(tot)); print("sample rows (fam,n,seed,l,#opt,#opt with RT edit,#opt with RR edit):", rows[:8])
