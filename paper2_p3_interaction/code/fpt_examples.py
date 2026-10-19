"""Small concrete examples quoted in notes/fpt_attempt.md."""
import itertools, random
from check_csp_view import pairs
from fpt_tests import *
def show(n, E, H):
    packed = {q for P in H for q in pairs(P)}
    o = opt(n, E); print("n", n, "E", sorted(map(sorted, E)), "H", H, "OPT", o, "|H|", len(H))
    for S in sol_sets(n, E, o):
        if len(S) == o: print("   opt S:", sorted(map(sorted, S)), "unpacked:", sorted(map(sorted, [q for q in S if q not in packed])))
# smallest instance where every optimum has an unpacked edit touching R
rng = random.Random(0); best = None
for n in range(4, 8):
    for _ in range(4000):
        E = rand_graph(rng, n); H = rand_packing(rng, n, E, True)
        if not H: continue
        packed = {q for P in H for q in pairs(P)}; Vh = {x for P in H for x in P}
        o = opt(n, E)
        if o <= len(H): continue
        sols = [S for S in sol_sets(n, E, o) if len(S) == o]
        if all(any(q not in packed and not q <= Vh for q in S) for S in sols):
            best = (n, E, H); break
    if best: break
show(*best)
