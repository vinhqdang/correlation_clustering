"""sparse chain search: l = 0, pair-maximal, I(H) path/tree, no reducible proper subfamily, max cluster size <= cmax in an optimum."""
import sys, random
from w6a_chain import *
def search2(seed, iters, m, r, cmax, maxn, path=True, ps=(0.05,0.1,0.2,0.3)):
    rng = random.Random(seed); found = []
    for it in range(iters):
        g = gen(m, r, rng.choice(ps), rng, path=path)
        if not g: continue
        n, adj, mem = g
        if n > maxn: continue
        if not all((adj[x] >> y & 1) and (adj[y] >> z & 1) and not (adj[x] >> z & 1) for x, y, z in mem): continue
        if not is_pair_maximal(adj, mem): continue
        o, part = opt(adj, True)
        if o != len(mem): continue
        if max(pc(c) for c in part) > cmax: continue
        T = list(range(len(mem))); W = members_vertices(mem, T)
        if not any(R_ok(adj, S, W) for S in tight_choices(adj, mem, T, W)): continue
        s = minred(adj, mem, T)
        if s == len(mem): found.append((n, sorted(edges_of(adj)), mem, [bin(c) for c in part]))
    return found
if __name__ == "__main__":
    seed, iters, m, r, cmax, maxn = map(int, sys.argv[1:7])
    f = search2(seed, iters, m, r, cmax, maxn)
    print("m", m, "r", r, "cmax", cmax, "found", len(f))
    for x in f[:3]: print(x)
