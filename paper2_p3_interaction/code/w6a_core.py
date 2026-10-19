"""Minimum joint-conflict core size among outside vertices for tight irreducible components (core = set Y of outside vertices such that no tight S_T satisfies all y in Y)."""
import sys, random, itertools, collections
from w6a_lib import *
from w6a_rule import *
from w6a_viol import y_ok
def main():
    seed, k, a, b = map(int, sys.argv[1:5]); rng = random.Random(seed); D = collections.Counter()
    for it in range(k):
        n = rng.randint(a, b); adj = random_instance(n, rng); H = maximal_pair_packing(adj, rng)
        if not H: continue
        for T in I_components(H):
            W = members_vertices(H, T)
            if opt(induced(adj, W)) > len(T): continue
            ch = tight_choices(adj, H, T, W)
            if any(R_ok(adj, S, W) for S in ch): continue
            outs = [y for y in range(n) if not W >> y & 1 and adj[y] & W]
            sat = {y: frozenset(i for i, S in enumerate(ch) if y_ok(adj, S, W, y)) for y in outs}
            core = None
            for r in range(1, len(outs) + 1):
                for Y in itertools.combinations(outs, r):
                    inter = set(range(len(ch)))
                    for y in Y: inter &= sat[y]
                    if not inter: core = r; break
                if core: break
            D[(len(T), core)] += 1
    print(sorted(D.items(), key=str))
main()
