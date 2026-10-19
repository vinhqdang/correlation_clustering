"""Search: forest I(H), maximal H, instances where EVERY optimum edits an unpacked pair (u,v), u,v in V(H),
whose P3s are at I(H)-distance >= 2 (or in different components).  Prints smallest found."""
import sys, random, itertools
import numpy as np
from w2a_lib import *
from w2a_bfk import forest_packing
def dist_ok(H, u, v):
    # u,v far: no P3 containing u adjacent (sharing a vertex) to a P3 containing v, and not in a common P3
    Pu = [i for i, P in enumerate(H) if u in P]; Pv = [i for i, P in enumerate(H) if v in P]
    for i in Pu:
        for j in Pv:
            if set(H[i]) & set(H[j]): return False
    return True
def run(trials, seed, nlo, nhi):
    rng = random.Random(seed); found = []; tested = 0
    for _ in range(trials):
        n = rng.randint(nlo, nhi); E = rand_graph(rng, n, rng.choice([.25, .35, .5]))
        H = forest_packing(rng, n, E)
        # make maximal in pair sense
        used = {q for P in H for q in pairs(P)}
        for P in induced_p3s(n, E):
            if not any(q in used for q in pairs(P)):
                H2 = H + [P]
                # keep forest property
                try:
                    from w2a_bfk import forest_packing as fp
                except: pass
                H = H2; used.update(pairs(P))
        if len(H) < 2: continue
        o = opt_cost(n, E); h = len(H)
        if o == h: continue
        Vh = sorted({x for P in H for x in P}); packed = used
        far = [(u, v) for u, v in itertools.combinations(Vh, 2) if fs(u, v) not in packed and dist_ok(H, u, v)]
        full = 1 << n; tested += 1
        N = [mask([y for y in range(n) if fs(x, y) in E]) for x in range(n)]
        al = np.ones(full, np.bool_)
        for S in range(full):
            ok = True
            for u, v in far:
                a, b = (S >> u) & 1, (S >> v) & 1
                if fs(u, v) in E:
                    if a != b: ok = False; break
                else:
                    if a and b: ok = False; break
            al[S] = ok
        if opt_cost(n, E, al) > o:
            found.append((n, o, h, sorted(map(sorted, E)), H)); 
    found.sort(key=lambda t: (t[0], len(t[3])))
    print("tested", tested, "found", len(found))
    for f in found[:2]: print(f)
if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
