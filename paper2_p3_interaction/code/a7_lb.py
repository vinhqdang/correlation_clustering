"""Audit 7: Corollary LB relaxation: min over soft-pair state vectors sigma, transitive on all-soft triples, of sum c_i + ceil(W(sigma)/2D).
Check (1) value <= OPT always; (2) is it exact at l=0 (value<=|H| => OPT==|H|)?"""
import sys, random, itertools, math
import numpy as np
import a7_lib as L
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2]); nr = [int(x) for x in sys.argv[3].split(',')]
tot = 0; lbviol = 0; nonexact = 0; ex = None; stronger = 0
for it in range(N):
    n = rng.choice(nr); A, H = L.gen_instance(n, rng)
    if not H or not L.is_pair_maximal(A, H): continue
    S = L.soft_of(n, H); Dl, Gm = L.delta_gamma(A, S)
    c = max(sum(1 for t in H if v in t) for v in range(n)); D = 2 * c
    sp = [p for p in L.prs(n) if S[p]]; m = len(sp); si = {p: k for k, p in enumerate(sp)}
    triples = [t for t in itertools.combinations(range(n), 3) if all(S[x, y] for x, y in itertools.combinations(t, 2))]
    opt = L.opt_subset(A); best = 10 ** 9
    for bits in range(1 << m):
        st = [(bits >> k) & 1 for k in range(m)]
        ok = True
        for t in triples:
            a, b, c3 = t
            s = st[si[(a, b)]] + st[si[(a, c3)]] + st[si[(b, c3)]]
            if s == 2: ok = False; break
        if not ok: continue
        ci = sum(1 for p in sp if st[si[p]] != int(A[p]))
        W = sum(len(Dl[p]) if st[si[p]] else len(Gm[p]) for p in sp)
        v = ci + math.ceil(W / (2 * D)); best = min(best, v)
    tot += 1
    if best > opt: lbviol += 1
    if best <= len(H) and opt > len(H): nonexact += 1; ex = (A.tolist(), H, best, opt)
    if best < opt: stronger += 1
print("instances", tot, "LB violations", lbviol, "l=0 non-exact (relax says |H|, OPT>|H|):", nonexact, "gap cases (value<OPT):", stronger)
print(ex)
