"""Test shrink rule: twin class X of R (H maximal) with |X| >= l+2 -> delete one vertex, same decision for budget |H|+l.
and: clique component (A_X empty) deletion."""
import sys, random, itertools
import numpy as np
from w2a_lib import *
def big_structured(rng, n):
    k = rng.randint(4, 6); E = rand_graph(rng, k, rng.choice([.3, .5, .7])); nxt = k; pool = list(range(k)); rng.shuffle(pool)
    while nxt < n:
        s = min(rng.randint(1, 5), n - nxt); X = list(range(nxt, nxt + s)); nxt += s
        A = [pool.pop() for _ in range(rng.randint(0, min(len(pool), 3)))]
        for x, y in itertools.combinations(X, 2): E.add(fs(x, y))
        for x in X:
            for z in A: E.add(fs(x, z))
    return E
def delete_vertex(n, E, H, x):
    mp = {v: (v if v < x else v - 1) for v in range(n) if v != x}
    E2 = {fs(mp[a], mp[b]) for a, b in map(tuple, E) if x not in (a, b)}
    H2 = [tuple(mp[v] for v in P) for P in H]
    return n - 1, E2, H2
def run(trials, seed):
    rng = random.Random(seed); st = dict(inst=0, shr=0, iso=0)
    for _ in range(trials):
        n = rng.randint(9, 11); E = big_structured(rng, n)
        if not induced_p3s(n, E): continue
        H = maximal_packing(n, E, rng, prefer={v for v in range(n) if v >= 5})
        Vh, R, cls = rclasses(n, E, H)
        if not cls: continue
        h = len(H); o = opt_cost(n, E); st['inst'] += 1
        for c in cls:
            x = c[-1]
            n2, E2, H2 = delete_vertex(n, E, H, x); o2 = opt_cost(n2, E2)
            A = [z for z in range(n) if fs(z, c[0]) in E and z not in c]
            for lp in range(0, 5):
                if len(c) >= lp + 2:
                    assert (o <= h + lp) == (o2 <= h + lp), ("shrink fails", n, sorted(map(sorted, E)), H, c, lp)
                    st['shr'] += 1
            if not A:
                assert o == o2 + 0 or o == o2, ("iso",)
                # delete whole class
                E3, n3, H3 = E, n, H
                for y in sorted(c, reverse=True): n3, E3, H3 = delete_vertex(n3, E3, H3, y)
                assert opt_cost(n3, E3) == o; st['iso'] += 1
    print(seed, st)
if __name__ == "__main__": run(int(sys.argv[1]), int(sys.argv[2]))
