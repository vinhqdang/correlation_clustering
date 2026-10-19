"""search chains: I(H) a path (or forest) of m members, l=0 (tight), pair-maximal, with NO reducible proper subfamily (all s(T)=|T|)."""
import sys, random, itertools, collections
from w6a_lib import *
from w6a_rule import *
from w6a_minred import minred

def gen(m, r, p, rng, path=True):
    mem = [(0, 1, 2)]; nv = 3
    for i in range(1, m):
        j = i - 1 if path else rng.randrange(i)
        # avoid sharing vertex that already shared twice? allow c<=2: pick vertex in mem[j] not already in 2 members
        cnt = collections.Counter(v for t in mem for v in t)
        cand = [v for v in mem[j] if cnt[v] < 2]
        if not cand: return None
        v = rng.choice(cand); a, b = nv, nv + 1; nv += 2
        t = [a, b]; t.insert(rng.randrange(3), v)
        mem.append(tuple(t))
    n = nv + r
    E = set(); forced = set()
    for (x, y, z) in mem:   # centre y
        E.add(key(x, y)); E.add(key(y, z)); forced |= {key(x, y), key(y, z), key(x, z)}
    # a tuple's middle element is centre: but vertex order was randomised -> ok, any order is a P3 (x,y,z)
    for a in range(n):
        for b in range(a + 1, n):
            if key(a, b) not in forced and rng.random() < p: E.add(key(a, b))
    return n, from_edges(n, E), mem

def search(seed, iters, m, r, maxn=12):
    rng = random.Random(seed); found = []
    for it in range(iters):
        g = gen(m, r, rng.choice([0.15, 0.3, 0.5, 0.7, 0.85]), rng, path=True)
        if not g: continue
        n, adj, mem = g
        if n > maxn: continue
        # H = mem (must be induced P3s with disjoint pairs: edges chosen so (x,y,z): xz non-edge is forced non-edge)
        ok = all((adj[x] >> y & 1) and (adj[y] >> z & 1) and not (adj[x] >> z & 1) for x, y, z in mem)
        if not ok or not is_pair_maximal(adj, mem): continue
        if opt(adj) != len(mem): continue
        T = list(range(len(mem)))
        W = members_vertices(mem, T)
        if not any(R_ok(adj, S, W) for S in tight_choices(adj, mem, T, W)): continue
        s = minred(adj, mem, T)
        if s == len(mem): found.append((n, sorted(edges_of(adj)), mem))
    return found

if __name__ == "__main__":
    seed, iters, m, r = map(int, sys.argv[1:5])
    f = search(seed, iters, m, r)
    print("m", m, "r", r, "found", len(f))
    for x in f[:3]: print(x)
