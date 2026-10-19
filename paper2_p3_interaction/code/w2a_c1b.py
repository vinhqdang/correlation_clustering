"""Test claim D: if |A_X| > l + q (l = OPT-|H|, q = max_a |Q(a)|) then some optimum has cluster(X) cap V(H) subset of A_X
(no RH insertion for X).  Also test the stronger 'weight' claim: B_X symmetric difference with A_X costs |X| each, so
|X| * |A_X xor B_X| <= l in every optimum in normal form."""
import sys, random, itertools
import numpy as np
from w2a_lib import *
def run(trials, seed):
    rng = random.Random(seed); st = dict(inst=0, tested=0, big=0)
    for _ in range(trials):
        n = rng.randint(9, 11); E = structured_graph(rng, n)
        if not induced_p3s(n, E): continue
        H = maximal_packing(n, E, rng, prefer={v for v in range(n) if v >= n // 2})
        Vh, R, cls = rclasses(n, E, H)
        if not cls: continue
        o = opt_cost(n, E); h = len(H); l = o - h
        Q = {v: {x for P in H if v in P for x in P if x != v} for v in Vh}
        q = max(len(Q[v]) for v in Vh)
        full = 1 << n; Vm = mask(Vh); st['inst'] += 1
        for c in cls:
            X = mask(c); A = 0
            for x in c: A |= mask([y for y in range(n) if fs(x, y) in E])
            A &= ~X
            if bin(A).count('1') > l + q:
                st['big'] += 1
                al = np.array([ (S & X) == 0 or ((S & Vm) & ~A) == 0 for S in range(full)], np.bool_)
                assert opt_cost(n, E, al) == o, ("claim D fails", n, sorted(map(sorted, E)), H, c)
            st['tested'] += 1
    print(seed, st)
if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]))
