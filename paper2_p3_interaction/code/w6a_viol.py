"""For tight irreducible components T: is there a forced violator y (every tight S_T violates y)?  Count components with NO forced violator
(every y individually satisfiable, but no S_T satisfies all y) -> 'joint conflict'."""
import sys, random, collections
from w6a_lib import *
from w6a_rule import *
def y_ok(adj, S, W, y):
    a2 = flip(adj, S); A = adj[y] & W
    if not A: return True
    return any(((a2[v] & W) | 1 << v) == A for v in range(len(adj)) if W >> v & 1)
def main():
    seed, k, a, b = map(int, sys.argv[1:5]); rng = random.Random(seed); st = collections.Counter(); ex = None
    for it in range(k):
        n = rng.randint(a, b); adj = random_instance(n, rng); H = maximal_pair_packing(adj, rng)
        if not H: continue
        for T in I_components(H):
            W = members_vertices(H, T)
            if opt(induced(adj, W)) > len(T): continue
            ch = tight_choices(adj, H, T, W)
            if any(R_ok(adj, S, W) for S in ch): continue
            st["tight_irr"] += 1
            outs = [y for y in range(n) if not W >> y & 1 and adj[y] & W]
            forced = [y for y in outs if not any(y_ok(adj, S, W, y) for S in ch)]
            if forced: st["has_forced"] += 1
            else:
                st["no_forced"] += 1
                if ex is None: ex = (n, sorted(edges_of(adj)), H, T, outs)
    print(dict(st)); print(ex)
if __name__ == "__main__": main()
