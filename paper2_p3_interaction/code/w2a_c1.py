"""Tests: twin lemma, C1 (no R-R edit), normal form, frozen classes, maximal-H facts. n=9..11."""
import sys, random, itertools
import numpy as np
from w2a_lib import *

def run(trials, seed, nlo=9, nhi=11):
    rng = random.Random(seed)
    st = dict(inst=0, gap=0, Rcls2=0, nontriv=0, twin=0, C1ok=0, frozen_tests=0, nf=0)
    for t in range(trials):
        n = rng.randint(nlo, nhi)
        E = structured_graph(rng, n) if rng.random() < .7 else rand_graph(rng, n, rng.choice([.2, .3, .45]))
        P0 = induced_p3s(n, E)
        if not P0: continue
        Rpref = set(range(n))  # prefer to use low-index (core) vertices: sort by index sum of high ones
        H = maximal_packing(n, E, rng, prefer={v for v in range(n) if v >= n // 2} if rng.random() < .6 else None)
        Vh, R, cls = rclasses(n, E, H)
        if len(R) < 2: continue
        # verify maximality + Lemma 4
        packed = {q for P in H for q in pairs(P)}
        for P in P0: assert any(q in packed for q in pairs(P))
        N = [mask([y for y in range(n) if fs(x, y) in E]) for x in range(n)]
        for r, s in itertools.combinations(R, 2):
            if fs(r, s) in E: assert N[r] | 1 << r == N[s] | 1 << s
            else: assert N[r] & N[s] == 0
        st['inst'] += 1
        o = opt_cost(n, E); h = len(H); st['gap'] += (o > h); st['Rcls2'] += (len(cls) >= 2)
        full = 1 << n
        # --- C1 normal form: each cluster meets R in empty or exactly one full class
        Rm = mask(R); clsm = {mask(c) for c in cls}
        allowed = np.array([ (S & Rm) == 0 or (S & Rm) in clsm for S in range(full)], np.bool_)
        o2 = opt_cost(n, E, allowed)
        assert o2 == o, ("C1 FAILS", n, sorted(map(sorted, E)), H)
        if o > h: st['nontriv'] += 1
        # --- twin lemma on R classes (adjacent true twins never separated in any optimum): separating them costs > OPT
        for c in cls:
            if len(c) >= 2:
                u, v = c[0], c[1]
                al = np.array([not (S >> u & 1 and S >> v & 1) for S in range(full)], np.bool_)
                assert opt_cost(n, E, al) > o, ("twin fails",)
                st['twin'] += 1
        # --- frozen: |X| > o-h  =>  optimum with cluster(X) = X u A_X exactly
        for c in cls:
            X = mask(c); A = 0
            for x in c: A |= N[x]
            A &= ~X
            if len(c) > o - h:
                al = np.array([ (S & X) == 0 or S == (X | A) for S in range(full)], np.bool_)
                assert opt_cost(n, E, al) == o, ("frozen fails", c)
                st['frozen_tests'] += 1
        # --- each V(H) vertex adjacent to <=1 class; non-edges inside A_X packed; edges a--A_X (a not in A_X u X) packed
        for a in Vh:
            assert sum(1 for c in cls if any(fs(a, x) in E for x in c)) <= 1
        for c in cls:
            A = {z for z in Vh if fs(z, c[0]) in E}
            for a, b in itertools.combinations(sorted(A), 2):
                if fs(a, b) not in E: assert fs(a, b) in packed
            for a in Vh - A:
                for z in A:
                    if fs(a, z) in E: assert fs(a, z) in packed
    print(seed, st)
if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]))
