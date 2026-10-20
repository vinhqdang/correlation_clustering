"""Audit 11: FAR_m repack (Prop 3.1) and Prop 5.2 (augmentation depth), independent re-derivation."""
from a11_lib import *
import itertools

def far(m):
    n = 2 * m + 1
    g = Gr(n, [(i, i + 1) for i in range(2 * m)] + [(1, 2 * m - 1)])
    H0 = [(2 * i, 2 * i + 1, 2 * i + 2) for i in range(m)]
    Hs = H0[:m - 1] + [(1, 2 * m - 1, 2 * m - 2)]
    return g, H0, Hs

for m in range(3, 10):
    g, H0, Hs = far(m)
    nu = nu_mod_ilp(g)
    opt = opt_ilp(g)[0]
    Is = I_graph(Hs)
    cyc = all(len(v) == 2 for v in Is.values())
    s = f'm={m} n={g.n} valid(H0)={is_packing(g,H0)} pm(H0)={pair_maximal(g,H0)} valid(H*)={is_packing(g,Hs)} pm(H*)={pair_maximal(g,Hs)} nu={nu} OPT={opt} I(H*)cycle={cyc} tw(I(H*))={tw_exact(Is)} tw(I(H0))={tw_exact(I_graph(H0))}'
    if g.n <= 11:
        o, opts = opt_and_optima(g)
        assert o == opt
        hs = [sorted(hard_of(Hs, S)) for S in opts]
        h0 = [sorted(hard_of(H0, S)) for S in opts]
        cov = {v for t in Hs for v in t}
        s += f' #opt={len(opts)} hard-edit counts wrt H*={sorted(set(len(h) for h in hs))} wrt H0={sorted(set(len(h) for h in h0))}'
        # all optima: is chord always deleted? is the hard edit wrt H* always at the uncovered vertex 2m?
        s += f' chord deleted in all={all((1,2*m-1) in S for S in opts)} all-H*-hard-edits-touch-uncovered={all(all((a not in cov) or (b not in cov) for a,b in h) for h in hs)}'
        # distance in I(H0) of ends of every hard edit
    print(s)

# Prop 5.2: path 0..2m+2, H_0 = {(2i+1,2i+2,2i+3)}; nu=m+1=OPT; unique max packing disjoint from H_0
print('--- Prop 5.2')
for m in range(2, 7):
    n = 2 * m + 3
    g = Gr(n, [(i, i + 1) for i in range(n - 1)])
    H0 = [(2 * i + 1, 2 * i + 2, 2 * i + 3) for i in range(m)]
    mx, mps = max_packings(g)
    opt = opt_ilp(g)[0]
    shared = [len(set(H0) & set(tuple(sorted((x[0], x[1], x[2]))) if False else x for x in mp)) for mp in mps]
    # members compared up to orientation
    norm = lambda t: (min(t[0], t[2]), t[1], max(t[0], t[2]))
    sh = [len({norm(t) for t in H0} & {norm(t) for t in mp}) for mp in mps]
    print(f'm={m}: pairmax(H0)={pair_maximal(g,H0)} valid={is_packing(g,H0)} nu={mx} OPT={opt} #max packings={len(mps)} common members with H0={sh} tw(I(H0))={tw_exact(I_graph(H0))}')
    # minimal r: smallest r such that some r-subset R of H0 and (r+1) pair-disjoint P3s avoiding pairs of H0-R exist
    P3 = g.p3()
    best = None
    for r in range(0, m + 1):
        found = False
        for R in itertools.combinations(range(m), r):
            keep = [H0[i] for i in range(m) if i not in R]
            used = {p for t in keep for p in pr(t)}
            cand = [t for t in P3 if not any(p in used for p in pr(t))]
            # need r+1 pair-disjoint among cand
            def rec(i, cur, u):
                if len(cur) == r + 1: return True
                if i == len(cand): return False
                if len(cur) + len(cand) - i < r + 1: return False
                t = cand[i]
                if not any(p in u for p in pr(t)):
                    for p in pr(t): u.add(p)
                    ok = rec(i + 1, cur + [t], u)
                    for p in pr(t): u.discard(p)
                    if ok: return True
                return rec(i + 1, cur, u)
            if rec(0, [], set()): found = True; break
        if found: best = r; break
    print('   minimal improving exchange size r =', best)
