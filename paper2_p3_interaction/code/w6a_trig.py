"""attach a C4 trigger (member M0=(p,q,r), s) to N(F) by a few random edges; want: pair-maximal, OPT=|H|+1, T (whole N(F)-component) irreducible, no reducible small subfamily."""
import sys, random, subprocess, itertools
from w6a_lib import *
from w6a_rule import *
from w6a_fam import build
def copt(adj):
    inp = str(len(adj)) + "\n" + " ".join(map(str, adj)) + "\n"
    return int(subprocess.run([__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "w6a_opt")], input=inp, capture_output=True, text=True).stdout)
def attach(adj0, n0, extra):
    n = n0 + 4; E = edges_of(adj0)
    p, q, r, s = n0, n0 + 1, n0 + 2, n0 + 3
    E |= {key(p, q), key(q, r), key(r, s), key(s, p)} | set(extra)
    return n, from_edges(n, E), (p, q, r)
def main():
    seed, iters, smax = map(int, sys.argv[1:4]); rng = random.Random(seed)
    tri = [(0, 1), (1, 2), (2, 0)]
    n0, adj0, H0, cl, sigma = build(tri)
    found = []
    for it in range(iters):
        k = rng.randint(1, 4)
        extra = set()
        for _ in range(k):
            t = rng.randrange(4); v = rng.randrange(n0); extra.add(key(n0 + t, v))
        n, adj, M0 = attach(adj0, n0, extra)
        H = H0 + [M0]
        if not all((adj[x] >> y & 1) and (adj[y] >> z & 1) and not (adj[x] >> z & 1) for x, y, z in H): continue
        if len(soft_pairs(H)) != 3 * len(H): continue
        if not is_pair_maximal(adj, H): continue
        comps = I_components(H)
        bad = False
        for T in comps:
            W = members_vertices(H, T)
            ch = tight_choices(adj, H, T, W)
            if len(T) > 1 and not ch: bad = True; break      # T must be tight
            if len(T) > 1 and any(R_ok(adj, S, W) for S in ch): bad = True; break   # T must be irreducible
        if bad: continue
        red = False
        for F in conn_subfamilies(H, smax):
            W = members_vertices(H, F)
            if pair_closed(H, F, W) and any(R_ok(adj, S, W) for S in tight_choices(adj, H, F, W)): red = True; break
        if red: continue
        o = copt(adj)
        found.append((o - len(H), sorted(extra)))
        if o - len(H) == 1:
            print("FOUND", sorted(extra), "n", n, flush=True); return
    print("none; stats", sorted(set(f[0] for f in found)), len(found))
if __name__ == "__main__": main()
