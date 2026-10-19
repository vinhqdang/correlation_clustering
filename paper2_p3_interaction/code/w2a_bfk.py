"""E2: c=1 consistency with van Bevern et al. (Rule 6.2 on single P3s, |H|<=2l if no rule applies).
   E6: Rule 6.2 on connected sub-packings (modules) of a FOREST I(H): soundness test."""
import sys, random, itertools
import numpy as np
from w2a_lib import *

def nbrs(n, E): return [{y for y in range(n) if fs(x, y) in E} for x in range(n)]
def parts_of(Wl):
    def rec(i, lab, m):
        if i == len(Wl): yield list(lab); return
        for l in range(m + 1):
            lab.append(l); yield from rec(i + 1, lab, max(m, l + 1)); lab.pop()
    yield from rec(0, [], 0)
def opt_solutions_of_W(Wl, E):
    """all optimal edit sets (frozenset of pairs) of G[W]"""
    best, sols = 10**9, []
    for lab in parts_of(Wl):
        S = frozenset(fs(a, b) for (i, a), (j, b) in itertools.combinations(list(enumerate(Wl)), 2)
                      if (lab[i] == lab[j]) != (fs(a, b) in E))
        if len(S) < best: best, sols = len(S), [(S, lab)]
        elif len(S) == best: sols.append((S, lab))
    return best, sols
def rule_cond(n, E, Wl, lab, N):
    W = set(Wl)
    for (i, u), (j, v) in itertools.combinations(list(enumerate(Wl)), 2):
        if lab[i] == lab[j]:
            if N[u] - W != N[v] - W: return False
        else:
            if (N[u] & N[v]) - W: return False
    return True

def E2(trials, seed):
    rng = random.Random(seed); st = dict(inst=0, rule_applied=0, noRule=0, noRule_pos=0)
    for _ in range(trials):
        n = rng.randint(8, 11); E = rand_graph(rng, n, rng.choice([.15, .25, .35, .5]))
        cand = induced_p3s(n, E); rng.shuffle(cand); used = set(); H = []
        for P in cand:
            if not (set(P) & used): H.append(P); used |= set(P)
        if not H: continue
        o = opt_cost(n, E); N = nbrs(n, E); st['inst'] += 1
        any_rule = False
        for P in H:
            Wl = list(P); tau, sols = opt_solutions_of_W(Wl, E)
            assert tau == 1
            for S, lab in sols:
                if rule_cond(n, E, Wl, lab, N):
                    any_rule = True; st['rule_applied'] += 1
                    E2g = set(E) ^ set(S)
                    assert opt_cost(n, E2g) == o - 1, ("Lemma 6.3 consequence fails", n, sorted(map(sorted, E)), P)
        if not any_rule:
            st['noRule'] += 1; l = o - len(H)
            assert 2 * l >= len(H), ("|H|<=2l fails", n, sorted(map(sorted, E)), H, o)
            st['noRule_pos'] += (len(H) > 0)
    print("E2", seed, st)

def forest_packing(rng, n, E):
    cand = induced_p3s(n, E); rng.shuffle(cand); used = set(); H = []
    def forest_ok(Hs):
        m = len(Hs); par = list(range(m))
        def f(x):
            while par[x] != x: par[x] = par[par[x]]; x = par[x]
            return x
        for i in range(m):
            for j in range(i):
                if set(Hs[i]) & set(Hs[j]):
                    a, b = f(i), f(j)
                    if a == b: return False
                    par[a] = b
        return True
    for P in cand:
        ps = pairs(P)
        if not any(q in used for q in ps) and forest_ok(H + [P]):
            H.append(P); used.update(ps)
    return H
def connected_subpackings(H, maxsize):
    m = len(H); adj = [[j for j in range(m) if j != i and set(H[i]) & set(H[j])] for i in range(m)]
    seen = set()
    for i in range(m):
        stack = [frozenset([i])]
        while stack:
            F = stack.pop()
            if F in seen: continue
            seen.add(F)
            if len(F) < maxsize:
                for a in F:
                    for b in adj[a]:
                        if b not in F: stack.append(F | {b})
    return seen
def E6(trials, seed):
    rng = random.Random(seed); st = dict(inst=0, mods=0, applied=0, nontriv_applied=0, excess_applied=0)
    for _ in range(trials):
        n = rng.randint(8, 11); E = rand_graph(rng, n, rng.choice([.2, .3, .4]))
        H = forest_packing(rng, n, E)
        if len(H) < 2: continue
        o = opt_cost(n, E); N = nbrs(n, E); st['inst'] += 1
        for F in connected_subpackings(H, 3):
            Wl = sorted({v for i in F for v in H[i]}); 
            if len(Wl) > 7: continue
            st['mods'] += 1
            tau, sols = opt_solutions_of_W(Wl, E)
            for S, lab in sols:
                if rule_cond(n, E, Wl, lab, N):
                    st['applied'] += 1; st['nontriv_applied'] += (len(F) >= 2); st['excess_applied'] += (tau > len(F))
                    E2g = set(E) ^ set(S)
                    assert opt_cost(n, E2g) == o - tau, ("module rule unsound", n, sorted(map(sorted, E)), H, sorted(F))
                    rest = [H[i] for i in range(len(H)) if i not in F]
                    for P in rest:  # remaining packing still induced P3s with the same pairs
                        assert not any(q in S for q in pairs(P))
                        a, b, c = P
                        assert fs(a, b) in E2g and fs(b, c) in E2g and fs(a, c) not in E2g
    print("E6", seed, st)
if __name__ == "__main__":
    E2(int(sys.argv[1]), int(sys.argv[2])); E6(int(sys.argv[3]), int(sys.argv[2]))
