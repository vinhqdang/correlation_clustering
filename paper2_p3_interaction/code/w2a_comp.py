"""E7: BFK Rule 6.2 / Lemma 6.5 applied to the COMPONENTS of I(H) (vertex-disjoint induced subgraphs G[V(C)]).
Checks (a) rule soundness OPT(G^S) = OPT - tau, (b) in a reduced instance, with k = OPT and l' = k - sum tau,
#components <= 2 l'  (BFK Lemma 3.2/6.5 at component level), and l' >= 0 always."""
import sys, random, itertools
import numpy as np
from w2a_lib import *
from w2a_bfk import nbrs, opt_solutions_of_W, rule_cond
def comps(H):
    m = len(H); par = list(range(m))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for i in range(m):
        for j in range(i):
            if set(H[i]) & set(H[j]): par[f(i)] = f(j)
    g = {}
    for i in range(m): g.setdefault(f(i), []).append(i)
    return list(g.values())
def cap_packing(rng, n, E, capv):
    cand = induced_p3s(n, E); rng.shuffle(cand); used = set(); H = []
    for P in cand:
        ps = pairs(P)
        if any(q in used for q in ps): continue
        H2 = H + [P]
        if all(len({v for i in c for v in H2[i]}) <= capv for c in comps(H2)):
            H = H2; used.update(ps)
    return H
def run(trials, seed):
    rng = random.Random(seed); st = dict(inst=0, applied=0, reduced=0, red_comp_pos=0, maxratio=0.0)
    for _ in range(trials):
        n = rng.randint(8, 11); E = rand_graph(rng, n, rng.choice([.2, .3, .4, .55]))
        H = cap_packing(rng, n, E, rng.choice([3, 5, 6, 7]))
        if not H: continue
        st['inst'] += 1
        cur = set(E); C = [sorted({v for i in c for v in H[i]}) for c in comps(H)]
        o = opt_cost(n, cur); tot_applied = 0
        changed = True
        while changed:
            changed = False
            N = nbrs(n, cur)
            for Wl in list(C):
                tau, sols = opt_solutions_of_W(Wl, cur)
                for S, lab in sols:
                    if rule_cond(n, cur, Wl, lab, N):
                        new = cur ^ set(S); assert opt_cost(n, new) == opt_cost(n, cur) - tau
                        cur = new; C.remove(Wl); tot_applied += tau; changed = True; st['applied'] += 1; break
                if changed: break
        k = opt_cost(n, cur)
        taus = [opt_solutions_of_W(Wl, cur)[0] for Wl in C]
        lp = k - sum(taus); assert lp >= 0
        st['reduced'] += 1
        if C:
            st['red_comp_pos'] += 1
            assert len(C) <= 2 * lp, ("Lemma 3.2 at component level fails", n, sorted(map(sorted, cur)), C, k, taus)
            st['maxratio'] = max(st['maxratio'], len(C) / max(lp, 1e-9))
    print(seed, st)
if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]))
