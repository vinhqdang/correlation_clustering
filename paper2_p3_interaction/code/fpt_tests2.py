"""T6: Lemma 5 (structure of R = V minus V(H) for maximal H) and
   T7: conjecture C1 (some optimum has no unpacked edit with both endpoints in R), n up to 8."""
import itertools, random
from check_csp_view import pairs, induced_p3s
from fpt_tests import fs, rand_graph, rand_packing, sol_sets, opt

def T6(trials, seed):
    rng = random.Random(seed); c = 0
    for _ in range(trials):
        n = rng.randint(4, 8); E = rand_graph(rng, n); H = rand_packing(rng, n, E, True)
        Vh = {x for P in H for x in P}; R = sorted(set(range(n)) - Vh)
        N = lambda v: {x for x in range(n) if fs(v, x) in E}
        for r, s in itertools.combinations(R, 2):
            if fs(r, s) in E: assert N(r) | {r} == N(s) | {s}
            else: assert not (N(r) & N(s))
            c += 1
    print("T6 ok on", c, "R-pairs")

def T7(trials, seed):
    rng = random.Random(seed); tested = bad = 0
    for _ in range(trials):
        n = rng.choice([7, 8]); E = rand_graph(rng, n); H = rand_packing(rng, n, E, True)
        if not H: continue
        Vh = {x for P in H for x in P}
        if n - len(Vh) < 2: continue
        packed = {q for P in H for q in pairs(P)}; o = opt(n, E)
        sols = [S for S in sol_sets(n, E, o) if len(S) == o]
        tested += 1
        if all(any(q not in packed and not (q & Vh) for q in S) for S in sols): bad += 1
    print("T7: tested", tested, "instances with |R|>=2; forced R-R edit in", bad)
if __name__ == "__main__":
    T6(600, 5); T7(500, 6)
